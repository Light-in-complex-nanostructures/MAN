%% ScriptQNM_web.m
% This code is used to calculate a QNM ; find its pole and then normalize the field
% It requires the use of COMSOL LiveLink.
% Beware denegerate states (and close poles)
% Beware the epsilon and mu of COMSOL material properties must depend on the frequency omega in rad/s
% At the end, the symmetry coefficient is only taken into account for the normalized field

clear all,close all
%   Matlab catches COMSOL exception during the execution of this code to
%   prevent program crash when simulation doesn't converge (which is expected
%   at some point), so to debug using COMSOL error messages, uncomment the
%   following line :
%debug_mode=1


%% %%%%%%%%%%%%%%%%% PARAMETERS TO BE DEFINED BY USER% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PHYSICAL PARAMETERS
% Initial guess pole wavelength to create the first 3 pulsations for iteration
lambda_guess=(920 - 50*1i)*1e-9;% in meter
% real(lambda_guess) is also used as omega_L for using COMSOL with a complex frequency, see appendix 1 in OE2013

%% NUMERICAL CRITICAL PARAMETERS to be controlled
save_file='data_rod_normalized.mat'; % .mat file name to save the data
save_model='save_model_rod_normalized.mph'; % .mph file to store normalized QNM
% intial frequencies
delta=0.001; % control the slight shift of the 3 initial guess frequencies omega=[omega_start, omega_start*(1+delta) omega_start*(1-delta)]
% Symmetry factor to account for symmetry used in calculation(has to be >1, =1 only if no symmetry at all)
sym_factor=4; % 1 (without sym) 2 (1 symmetry plane) or 4 (2 symmetry planes)

%% COMSOL PARAMETERS
% COMSOL model file name
model_name='models/nanorod_OE2013.mph';% if you use an other model sheet, check line 65
% Name of the frequency parameter used in the tag "materials" in COMSOL sheet (epsilon(omega) and mu(omega)) in rad/s
omega_var_name='omega';
% Material of the background (the background has to be homogeneous) in the tag "material 1" of COMSOL sheet
background_material='mat1';
% maximum number of COMSOL-Solver iterations
max_solv_iter=100; % Since the matrix is singular, close to the pole, COMSOL may perform a huge amount of iterations for nothing

%% NUMERICAL PARAMETERS
% Maximum number of iterations to iteratively compute the pole (a typical value is 5, above 10 the research is not converging)
QNM_ite=10;
% Coordinate of the point at which the field is calculated for the iterative pole search ((not the dipole location)
x_sample=0; y_sample=0; z_sample=-80e-9; % (80 nm) in the distance-unit used in COMSOL model sheet
% component of the fields (Ex, Ey, Ez, Hx, Hy, Hz) used to find the pole (Ez in OExpress 2013)
tested_field_comp='Ez';
% window & resolution to visualize the field on xOz plane
window=[0 100 -150 150]*1e-9;% [xmin xmax zmin zmax] (m)
resolution=0.5e-9;% (0.5nm) pixel size
%%%%%%%%%%%%%%%%%%% PARAMETERS TO BE DEFINED BY USER (finished) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% %%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZE POLE-SEARCH LOOP
import com.comsol.model.* 		% Load APIs
import com.comsol.model.util.*	% idem
ModelUtil.showProgress(true) ;  % Progress bar
% COMSOL model tags used in the script; (but be aware of them in case of a bug)
% DO NOT CHANGE, except if you have several geometries or studies ...
total_field_study_name='std1';bg_field_study_name='std2';geometry_name='geom1';mesh_name='mesh1';electric_point_dipole_name='epd1';solver_name='sol1';
% evaluation point
eval_point=[x_sample;y_sample;z_sample];
% Create the initial triplets to estimate the pole
c=299792458;% Speed of light
omega_start=2*pi*299792458/lambda_guess; % Initial guess value
omega_set=omega_start;
tested_field_set=0;
omega=omega_start;
pole_estimate=omega_start;
% Preparation of figures
h1=figure;
% Window for normalized field visualisation (xOz plane)
xp=window(1):resolution:window(2);% create the x axis 
zp=window(3):resolution:window(4);% create the z axis
[Xp,Zp]=meshgrid(xp,zp);
XYZ=[Xp(:),0*Xp(:),Zp(:)];% create coordinates to which we evaluate the field map
% Create a custom colormap
v=linspace(0,1,65);v=v(2:end);iv=fliplr(v);
z=zeros(1,64);
o=ones(1,64);
r=[ z, z, v, o];
g=[iv, z, z, v];
b=[ o,iv, z, z];
cmap=[r',g',b'];
clear p x_sample y_sample z_sample Xp Zp v iv z o r g b

%% POLE-SEARCH LOOP
for iter=1:QNM_ite
    tic ;
    fprintf('Iteration %d ...\n ',iter);
    %% Load COMSOL model and adapt to complex frequencies
    model = mphload(model_name); %we reload the model to refresh it
    setCOMSOL_ComplexFreq(model,omega_var_name,real(omega));
    model.result.numerical.create('global_eval_QNM', 'Global');%tool to evaluate COMSOL quantities
    % limit the number of solver iterations  (only possible if the solver 'sol1' is already
    % defined in the model sheet, and is the only solver defined)
    model.sol(solver_name).feature('s1').feature('i1').set('maxlinit', max_solv_iter);
    
    %% Run the COMSOL simulation at 'omega'
    fprintf('\tCOMSOL field calculations on pole frequency estimate\n');
    % set the frequency of COMSOL computation
    model.param.set(omega_var_name,num2str(omega,'%12.15e'));
    model.geom(geometry_name).run();
    model.mesh(mesh_name).run();
    
    % Run the simulation
    if ~exist('debug_mode','var') % Matlab catches COMSOL exception to avoid program crash
        try
            model.study(total_field_study_name).run(); % at this stage,the Maxwell eqs have been solved
        catch
            fprintf('\n Iteration %s did not converge on users criteria: QNM pole reached (possibly)\n',num2str(iter));
            fprintf('\n\t Pole : %1.15e + %1.15e I \n',real(omega), imag(omega));
            toc
            break
        end
    else
        model.study(total_field_study_name).run(); % at this stage,the Maxwell eqs have been solved
    end
    
    % get the field on the point at which the field is calculated for the iterative pole search
    [temp]=mphinterp(model,{['emw.' tested_field_comp]},'coord',eval_point,'Complexout','on');
    tested_field_tot(iter)=temp; clear t_field
    clear temp
    
    % Find information on the electric-dipole in COMSOL model sheet
    % We check that user defined the point dipole by a real dipole moment (and not by amplitude & direction).
    if strcmp( model.physics('emw').feature(electric_point_dipole_name).getString('DipoleSpecification'),'DipoleMoment')
        dip_mom_str=model.physics('emw').feature(electric_point_dipole_name).getStringArray('pI');
        model.result.numerical('global_eval_QNM').set('expr',dip_mom_str(3));% Jz
        dipole_mom(3)=model.result.numerical('global_eval_QNM').getReal();
        model.result.numerical('global_eval_QNM').set('expr',dip_mom_str(2));% Jy
        dipole_mom(2)=model.result.numerical('global_eval_QNM').getReal();
        model.result.numerical('global_eval_QNM').set('expr',dip_mom_str(1));% Jz
        dipole_mom(1)=model.result.numerical('global_eval_QNM').getReal();
        %We retrieve electric point dipole position
        num_dipole_entity= model.physics('emw').feature(electric_point_dipole_name).selection().entities(0);%COMSOL point entity the dipole is attached to
        dipole_pos= mphgetcoords(model,geometry_name,'point',num_dipole_entity);% position of this point
    else
        fprintf('\n Error : In COMSOL model sheet, electric dipole must be defined by a vector only\n');
        break
    end
    
    % get the field on the dipole
    ddd = mpheval(model,{'Ex','Ey','Ez'},'edim','point','selection',num_dipole_entity,'Complexout','on');
    E_source(iter,:)=[ddd.d1(:),ddd.d2(:),ddd.d3(:)];
    clear dip_mom_str
    
    % get the field on the xOz plane (on each point of the future field map)
    [total_field]=mphinterp(model,{['emw.' tested_field_comp]},'coord',XYZ.','Complexout','on');
    
    %% Run the COMSOL simulation at 'omega' without the resonator to obtain the background field
    fprintf('\tCOMSOL background calculations on pole frequency estimate\n');
    % We replace nanoresonator's materials by background's
    Permittivity_bg=char(model.material(background_material).propertyGroup('def').getString('relpermittivity'));
    Permeability_bg=char(model.material(background_material).propertyGroup('def').getString('relpermeability'));
    Nanores_materials=model.material.tags();%find the ID of all materials
    for i=1:numel(Nanores_materials) % change permeability and permittivity of all materials
        Permittivity{i}=char(model.material(Nanores_materials(i)).propertyGroup('def').getString('relpermittivity'));
        Permeability{i}=char(model.material(Nanores_materials(i)).propertyGroup('def').getString('relpermeability'));
        if ~strcmp(Nanores_materials(i),background_material)
            model.material(Nanores_materials(i)).propertyGroup('def').set('relpermittivity', {Permittivity_bg});
            model.material(Nanores_materials(i)).propertyGroup('def').set('relpermeability', {Permeability_bg});
        end
    end
    
    % We create a 2nd study to store data for the background field calculations
    model.study.create(bg_field_study_name);
    model.study(bg_field_study_name).create('freq', 'Frequency');
    model.study(bg_field_study_name).feature('freq').activate('emw', true);
    model.study(bg_field_study_name).feature('freq').set('plist', 'frequency_L_QNM');
    model.study(bg_field_study_name).feature('freq').set('punit', 'Hz');
    
    % Run the simulation
    model.study(bg_field_study_name).run(); % at this stage,the Maxwell eqs have been solved
    
    % Reset the resonator materials
    for i=1:numel(Nanores_materials) % change permeability and permittivity of all materials
        if ~strcmp(Nanores_materials(i),background_material)
            model.material(Nanores_materials(i)).propertyGroup('def').set('relpermittivity', {Permittivity{i}});
            model.material(Nanores_materials(i)).propertyGroup('def').set('relpermeability', {Permeability{i}});
        end
    end
    
    % get the field on the evaluation point
    [temp]=mphinterp(model,{['emw.' tested_field_comp]},'dataset','dset2','coord',eval_point,'Complexout','on');
    tested_field_bg(iter)=temp;
    clear temp
    
    % get the E field on the dipole
    ddd = mpheval(model,{'Ex','Ey','Ez'},'dataset','dset2','edim','point','selection',num_dipole_entity,'Complexout','on');
    E_source_bg(iter,:)=[ddd.d1(:),ddd.d2(:),ddd.d3(:)];
    
    % get the field on the xOz plane (on each point of the future field map)
    [background_field]=mphinterp(model,{['emw.' tested_field_comp]},'dataset','dset2','coord',XYZ.','Complexout','on');
    
    clear Ex Ey Ez Hx Hy Hz t_field
    
    %% Compute new set of (3) interpolating frequencies
    % generate the new frequency from triplet (see OE by Qiang Bai)
    tested_field_set(end)=tested_field_tot(iter);% we add the newly calculated field in the set
    [omega_set,tested_field_set]=omega_generation(omega_set,tested_field_set,delta);
    omega=omega_set(end);% the last value of the new omega set is the next to be calculated (new pole estimation)
    pole_estimate(iter+1)=omega;  % save the frequency at each calculation
    time_cal(iter)=toc;
    
    %% Plot field value for each iteration (|Field| must diverge when we reach the pole)
    figure(h1);
    if iter>3
        % PLOT the field divergence
        subplot(1,2,1)
        semilogy(1:iter,1./abs(tested_field_tot),'-mo',1:3,1./abs(tested_field_tot(1:3)),'bo-');
        xlabel('Iteration number');ylabel(['1/|' tested_field_comp '| at evaluation point']);title('Field divergence as one approaches the pole');
        ax=gca;set(ax,'XTick',1:iter);clear ax
        legend('intermediate estimates','initial triplet','Location','southwest');
        % PLOT the pole frequency estimations convergence
        subplot(2,2,2)
        plot(real(pole_estimate),'mo-'); hold on; plot(real(pole_estimate(1:3)),'bo-'); hold off ;
        ax=gca;set(ax,'XTick',1:iter);clear ax
        title('Pole real part (rad/s)');
        subplot(2,2,4)
        plot(imag(pole_estimate),'mo-'); hold on; plot(imag(pole_estimate(1:3)),'bo-'); hold off ;
        ax=gca;set(ax,'XTick',1:iter);clear ax
        xlabel('Iteration number');title('Pole imaginary part (rad/s)');
        drawnow
    else
        % PLOT the field divergence
        subplot(1,2,1)
        semilogy(1:iter,1./abs(tested_field_tot),'bo-');
        xlabel('Iteration number');ylabel(['1/|' tested_field_comp '| at evaluation point']);title('Field divergence as one approaches the pole');
        ax=gca;set(ax,'XTick',1:iter);clear ax
        %legend('initial triplet','Location','southwest');
        % PLOT the pole frequency estimations convergence
        subplot(2,2,2)
        plot(real(pole_estimate),'bo-');
        ax=gca;set(ax,'XTick',1:iter);clear ax
        title('Pole real part (rad/s)');
        %legend('intermediate estimates','initial triplet','pole','Location','best');
        subplot(2,2,4)
        plot(imag(pole_estimate),'bo-');
        ax=gca;set(ax,'XTick',1:iter);clear ax
        xlabel('Iteration number');title('Pole imaginary part (rad/s)');
        drawnow
    end
    
    %% Save data
    clear Nanores_materials Permeability_bg Permittivity_bg Permeability Permittivity ddd
    save(save_file,'-regexp','^(?!(model|ans|h1)$).');
    mphsave(model,save_model);% We save the model with all the data
end
if iter<4 % If the simulation didn't even pass the initial triplet => stop all
    return
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(h1);
subplot(2,2,2);hold on;plot(iter,real(pole_estimate(end)),'ro','LineWidth',2);
subplot(2,2,4);hold on;plot(iter,imag(pole_estimate(end)),'ro','LineWidth',2);

%% Normalization coefficients
% get the scattered field on the dipole
% Since it seems we cannot use 'mpheval' with a 'join' dataset on COMSOL, we have to
% use previously evualed total and background field
E_source_scatt=E_source-E_source_bg;
% calculate the normalization coefficient (cf OE 2013)
normal_coeff=sqrt(1i*(omega-pole_estimate(1:(end-1)))./(dot(repmat(dipole_mom,iter-1,1),E_source_scatt,2)).');

%% PLOT the normalized field (tested component) on evaluation point to show the convergence of normalization process
tested_field_scatt=tested_field_tot-tested_field_bg;
tested_field_normalized=normal_coeff.*tested_field_scatt./sqrt(sym_factor);
normal_coeff=normal_coeff.*sign(real(tested_field_normalized));
tested_field_normalized=tested_field_normalized.*sign(real(tested_field_normalized));
figure;
subplot(2,1,1)
plot(real(tested_field_normalized),'mo-'); hold on; plot(real(tested_field_normalized(1:3)),'bo-');
xlabel('Iteration number');ylabel(['Re($\tilde{' tested_field_comp '}$)'],'Interpreter','LaTex');
title(['Re($\tilde{' tested_field_comp '}$) at evaluation point'],'Interpreter','LaTex');
ax=gca;set(ax,'XTick',1:iter);clear ax
legend('intermediate estimates','initial triplet','Location','best');
subplot(2,1,2)
plot(imag(tested_field_normalized),'mo-'); hold on; plot(imag(tested_field_normalized(1:3)),'bo-');
xlabel('Iteration number');ylabel(['Im($\tilde{' tested_field_comp '}$)'],'Interpreter','LaTex');
title(['Im($\tilde{' tested_field_comp '}$) at evaluation point'],'Interpreter','LaTex');
ax=gca;set(ax,'XTick',1:iter);clear ax
fprintf('\n Normalized value of the QNM at evaluation point\n');
fprintf('\n\t value : %1.15e + %1.15e I \n',real(tested_field_normalized(end)), imag(tested_field_normalized(end)));

%% PLOT the normalized field (tested component) on the plane xOz
% we reshape the map from a point list to an image
total_field=reshape(total_field,length(zp),length(xp));
background_field=reshape(background_field,length(zp),length(xp));
scatt_field=total_field-background_field;
normalized_field=normal_coeff(end)*scatt_field./sqrt(sym_factor);% Normalize the field
clear XYZ
figure
subplot(1,2,1);
imagesc(xp*1e9,zp*1e9,real(normalized_field)); colormap(cmap);colorbar;axis image; axis xy; caxis([-1 1]*12e15);
hold on
plot([0 15e-9 15e-9 0],[-50e-9 -50e-9 50e-9 50e-9],':w','Linewidth',2)
xlabel('x (nm)');ylabel('z (nm)');title(['Re($\tilde{' tested_field_comp '}$) on plane xOz'],'Interpreter','LaTex');
subplot(1,2,2);
imagesc(xp*1e9,zp*1e9,imag(normalized_field)); colormap(cmap);colorbar;axis image; axis xy; caxis([-1 1]*2.5e14);
hold on
plot([0 15e-9 15e-9 0],[-50e-9 -50e-9 50e-9 50e-9],'r','Linewidth',2)
xlabel('x (nm)');ylabel('z (nm)');title(['Im($\tilde{' tested_field_comp '}$) on plane xOz'],'Interpreter','LaTex');

%% SAVE DATA
model=mphload(save_model);
model.param.set('omega_N',num2str(pole_estimate(end-1),'%12.15e'));% send normalisation frequency back to save_model
model.param.set('pole',num2str(pole_estimate(end),'%12.15e'));% send QNM pole frequency back 
model.param.set('N_coef',num2str(normal_coeff(end),'%12.15e'));% send normalisation coefficient back to save_model 
model.param.set('sym_factor',num2str(sym_factor));% send symmetry factor back to save_model (not taken into account in N_coef !)
mphsave(model,save_model);
clear h1 model ans num_dipole_entity Nanores_materials Permeability_bg Permittivity_bg ddd  ...
    to_plot_cp_bg  i j coeff_str h1 
save(save_file);

