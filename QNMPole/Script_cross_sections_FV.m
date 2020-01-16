%% Script_cross_sections_FV.m
% This code is used to calculate the scattering and absorption cross_sections
% of a nanoresonator with full-vectorial calculations, using COMSOL.
% It requires the use of COMSOL LiveLink.

%% TO READ
% -This script calculates the cross-sections for a plane wave illumination PW+k= Ez.exp(i(+k.x-w.t) defined in line 85
% -This program cannot exploit symmetries (but can easily be modified to do it)
% -There mus not be any sources in the model sheet (the script automatically provides the model with the plane wave to perform the calculations)
% -Be carreful; this program has initially been develop to compare with QNM
% calculations, there might be some reminiscent in the way the script is written

% clear all,close all



%% %%%%%%%%%%%%%%%%% PARAMETERS TO BE DEFINED BY USER% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NUMERICAL PARAMETERS
% Save data
save_file='data_cross_sections_FV.mat'; % to save data
% Wavelength span for cross_section calculations
wavelength=linspace(700,1200,50);% nm

%% COMSOL PARAMETERS
% Name of the model (without symmetries exploited)
model_name='models/nanorod_OE2013_full.mph';
% Name of the frequency parameter used in the tag "materials" in COMSOL sheet (epsilon(omega) and mu(omega)). in rad/s
omega_var_name='omega';
% Material of the background (the background has to be homogeneous) in the tag "material 1" of COMSOL sheet)
background_material='mat1';
% Material of the resonator (the resonator has to be homogeneous) in the tag "material 2" of COMSOL sheet)
resonator_material='mat2';
% Explicit selection of boundaries (in COMSOL model sheet) on which we will
% perform Poyting flux integration to calculate scattering cross-sections
% (this allows the program to find the flux integration surface in the list
% of all existing geometrical boundaries in the model sheet)
int_surface_selection='sel1';





%% %%%%%%%%%%%%%%%%% CROSS-SECTION FULLY VECTORIAL CALCULATION WITH COMSOL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
fprintf('\nInitialisation COMSOL simulations ...\n');
%% INITIALISATION 
% some constant for computation
c=2.99792458e8;
eps0=8.854187817e-12;
mu0=4*pi*1e-7;
import com.comsol.model.* 		% Load APIs
import com.comsol.model.util.*	% idem
ModelUtil.showProgress(true) ;  % Progress bar
% COMSOL model tags used in the script; (but be aware of them in case of a bug)
% DO NOT CHANGE, except if you have several geometries or studies ...
study_name='std1';geometry_name='geom1';mesh_name='mesh1';electric_point_dipole_name='epd1';solver_name='sol1';

% load the model sheet
model = mphload(model_name);
omegas = 2*pi*c./wavelength*1e9;

% Retrieve the domains numbers composing the resonnator (to integrate onto them)
doms=1:model.geom(geometry_name).getNDomains();% list of geometrical domains in the model sheet
doms_bg=model.material(background_material).selection.entities(3);% list of geometrical domains composing the background
doms_res=setdiff(doms,doms_bg);% list of geometrical domains of the resonnator (domains that are not background)

% Retrieve the permitivity (expressions) for all resonator and background
% material
eps_bg_str=char(model.material(background_material).propertyGroup('def').getString('relpermittivity'));
eps_res_str=char(model.material(resonator_material).propertyGroup('def').getString('relpermittivity'));
model.study(study_name).feature('freq').set('plist', [omega_var_name '/2/pi']);
% Set COMSOL to scattered field formulation
model.physics('emw').prop('BackgroundField').set('SolveFor', 'scatteredField');



%% CROSS-SECTIONS CALCULATIONS
h=waitbar(0,'Calculate cross-sections with COMSOL ...');
fprintf('Calculating cross-sections with COMSOL ...\n');
for plo=1:length(wavelength)
    tic
    % Set the frequency in COMSOL model sheet
    model.param.set(omega_var_name,num2str(omegas(plo),'%12.15e')); % make sure 'omega' is defined in the 'parameter' of model sheet
    
    % Set the plane wave expression
    model.physics('emw').prop('BackgroundField').set('Eb',...
        {'0' '0' ['exp(i*sqrt(' eps_bg_str ')*' omega_var_name '/c_const*x)']});% Ez*exp(+i*k*x)
   
    % Calulate Poyting flux of the plane wave (amplitude 1 [V/m])
    model.param.set('PF_pw',['real(sqrt(' num2str(eps0/mu0,'%12.15e') '*' eps_bg_str '))/2']);   
    
    % Run the simulation
    model.geom(geometry_name).run();
    model.mesh(mesh_name).run();
    model.study(study_name).run();
    
    % Calculate the absorption cross section
    sigma_abs(plo)=omegas(plo)*eps0*mphint2(model,['-1/2/PF_pw*(emw.normE^2)*imag(' eps_res_str ')'],3,'selection',doms_res);% plane wave amplitude is 1[V/m]
    % Calculate the scattering cross section
    sigma_scat(plo)=mphint2(model,'(emw.relPoavx*nx+emw.relPoavy*ny+emw.relPoavz*nz)/PF_pw',2,'selection',[geometry_name '_' int_surface_selection]);
    
    waitbar(plo/length(wavelength),h,'Calculate cross-sections with COMSOL ...');
    toc
end
close(h)

%% Plot the cross-sections
h1=figure;
subplot(1,2,1)
plot(wavelength,sigma_abs,'m--');
title('Absorption cross-section');
xlabel('lambda (nm)');ylabel('$\sigma_A (m^2)$','Interpreter','LaTex');
subplot(1,2,2)
plot(wavelength,sigma_scat,'m--');
title('Scattering cross-section');
xlabel('lambda (nm)');ylabel('$\sigma_S (m^2)$','Interpreter','LaTex');

clear h i ans model h1 eps_bg_str eps_res_str omegas choices
save(save_file)