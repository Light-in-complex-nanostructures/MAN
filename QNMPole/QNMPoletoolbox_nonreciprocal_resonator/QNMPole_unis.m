clear
%% USER DEFINED PARAMETERS 1 
% Initial guess pole wavelength to create the first 3 pulsations for iteration
% real(lambda_guess) is also used as omega_L for using COMSOL with a complex frequency, see appendix 1 in OE2013
lambda_guess=3e8/((8.65+0.055i)*1e9);% in meter
cms.delta=0.005; % control the slight shift of the 3 initial guess frequencies omega=[omega_start, omega_start*(1+delta) omega_start*(1-delta)]
% symmetry factor to account for symmetry used in calculation(has to be >1, =1 only if no symmetry at all)
cms.sym_factor=1; % 1 (without sym) 2 (1 symmetry plane) or 4 (2 symmetry planes)
% Input COMSOL model file name
cms.model_name='...\QNMPole_wire_YIG.mph';
% Files to save the computed results 
cms.save_file=  '...\mode.mat';
cms.save_model= '...\mode.mph';
% Name of the frequency parameter used in the tag "materials" in COMSOL sheet (epsilon(omega) and mu(omega)) in rad/s
cms.omega_var_name='omega';
% Material of the background (the background has to be homogeneous) in the tag "material 1" of COMSOL sheet
cms.background_material='wee2';

%% USER DEFINED PARAMETERS 2: ITERATION CONTROL
% Maximum number of iterations to iteratively compute the pole (a typical value is 5, above 10 the research is not converging)
cms.QNM_ite=15;
% Coordinate of the point at which the field is calculated for the iterative pole search ((not the dipole location)
cms.x_sample=-1.3e-2; cms.y_sample=5.2e-3;  % (80 nm) in the distance-unit used in COMSOL model sheet
% component of the fields (Ex, Ey, Ez, Hx, Hy, Hz) used to find the pole (Ez in OExpress 2013)
cms.tested_field_comp='Ez';

%% USER DEFINED PARAMETERS 3: ADVANCED
import com.comsol.model.* 		% Load APIs
import com.comsol.model.util.*	% idem
ModelUtil.showProgress(true) ;  % Progress bar
% COMSOL model tags used in the script; (but be aware of them in case of a bug)
% DO NOT CHANGE, except if you have several geometries or studies ...
cms.total_field_study_name='std1';cms.bg_field_study_name='std2';cms.geometry_name='geom1';cms.mesh_name='mesh1';cms.electric_point_dipole_name='lco1';cms.solver_name='sol1';
cms.res_field_study_name='std3';
% evaluation point
cms.eval_point=[cms.x_sample;cms.y_sample];
% Create the initial triplets to estimate the pole
c=299792458;% Speed of light
LH.omega_start=2*pi*c/lambda_guess; % Initial guess value
LH.omega_set=LH.omega_start;
LH.tested_field_set=0;
LH.omega=LH.omega_start;
LH.pole_estimate=LH.omega_start;
LH.h1=figure;

%% POLE-SEARCH LOOP
LH.E_source=zeros(1,3); LH.time_cal=0; LH.tested_field_tot=0;
LH.itercrash=0; 
for iter=1:cms.QNM_ite
    tic ;
    fprintf('Iteration %d ...\n ',iter);
    % Load COMSOL model and adapt to complex frequencies
    model = mphload(cms.model_name); % We reload the model to refresh it
    [LH.E_source, LH.tested_field_set, LH.tested_field_tot, LH.pole_estimate, LH.time_cal, LH.itercrash, LH.omega_set, LH.omega]=iteration( ...
        model, LH.omega, cms.omega_var_name, cms.geometry_name, cms.mesh_name, LH.tested_field_tot, LH.h1,...
        cms.total_field_study_name, cms.electric_point_dipole_name, cms.tested_field_comp, cms.eval_point, iter, ...
        LH.E_source, LH.tested_field_set, LH.pole_estimate, LH.time_cal, LH.omega_set, cms.delta, cms.save_model, cms.QNM_ite);
    if (LH.itercrash==1) % Check whether comsol has crashed
        break
    end
    % Save data
    clear Nanores_materials Permeability_bg Permittivity_bg Permeability Permittivity ddd
    save(cms.save_file,'-regexp','^(?!(model|ans|h1)$).');
end
if iter<4 % If the simulation didn't even pass the initial triplet => stop all
    return
end
%% POST-PROCESSING
% compute the background field
[cms.E_source_bg, cms.tested_field_bg]=bgcompute(cms, LH, iter);

%% NORMALIZATION COEFFICIENTS
% compute the normalization factor
[cms, LH]= QNMnorm(cms, LH, iter);

%% COMPUTE THE RESULTS FOR THE RECIPROCAL SYSTEM
[RH]=recip(cms, LH, iter);
[cms, RH]= QNMnorm(cms, RH, iter);
[signQNM]= QNMsign(cms, LH, RH, iter);
%% SAVE THE RESULTS
writemph(cms, LH, RH, signQNM)
clear model
save(cms.save_file)
%% 
% real(data1(emw.Ez*N_coef)*data2(emw.Ez*N_coefR))
