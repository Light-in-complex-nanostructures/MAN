%% Script_cross_sections_QNM.m
% This Matlab script calculates the scattering and absorption cross sections of a nanoresonator 
% for a plane wave illumination propagating along x and polarized along z
% Before, the QNM must have been calculated with 'ScriptQNM_web.m'
% The script relies on the QNM formalism published in Bai et al., Opt. Express 21, 27371-82 (2013).
% It requires the use of COMSOL LiveLink.

%% TO READ : BEWARE SYMMETRIES !
% The symmetry factor is defined as 1 (without sym) 2 (1 symmetry plane) or 4 (2 symmetry planes), see 'ScriptQNM_web.m' 
% If the symmetry factor is >=2, the program will assume there is the symmetry yOz plane (default: the nanorod model used in the demonstration code)
% NB1:  With symmetry, the overlap integral can thus only be calculated in one half (quarter) of the resonator
% NB2:  To calculate the integral over the entire resonator, a naive idea would be to say : because of symmetry, 
%       sub-integral onto each half (quarter) of the resonnator are all identical. But the plane wave is not symmetric by the yOz plane.
%       Thus we have to perform one integral for the right half of the resonator with plane wave +k.x, and simulate the integral on the 
%       left half of the QNM with PW +k.x by integrating the right half of the QNM with PW -k.x;
%       I[x=-a:+a](QNM(x).PW+k(x))=I[-a,0](QNM(x).PW+k(x))+I[0,+a](QNM(x).PW+k(x))
%                       =I[0,a](QNM(-x).PW+k(-x))+I[0,+a](QNM(x).PW+k(x))
%       if symmetry : QNM(-x)=QNM(x). 
%       Since PW+k(x)=Ez.exp(ik.x),PW-k(x)=Ez.exp(-ik.x)=PW+k(-x)
%       So  I[x=-a:+a](QNM(x).PW(x))= I[x=0:+a](QNM(x).PW(x))+ I[x=0:+a](QNM(x).PW-k(x))
clear all,close all


%% %%%%%%%%%%%%%%%%% PARAMETERS TO BE DEFINED BY USER% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NUMERICAL PARAMETERS
% Save data
save_file='data_cross_sections_QNM.mat'; % to save data
%
% Data from QNM calculation
normalisation_file='data_rod_normalized.mat';% all the other parameters needed will be loaded from this file (if the file hasn't been altered !)
% Material of the resonator (the resonator has to be homogeneous) in the tag "material 2" of COMSOL sheet)
resonator_material='mat2';
% Wavelength span for cross_section calculations
wavelength=700:5:1200;% nm





%% %%%%%%%%%%%%%%%%% CROSS-SECTION WITH QNM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALISATION
fprintf('\nInitialisation overlap integrals ...\n');
% some constant for computation
c=2.99792458e8;                 % speed of light in vacuum
eps0=8.854187817e-12;           % permittivity of vacuum
mu0=4*pi*1e-7;                  % permeability of vacuum
import com.comsol.model.* 		% Load APIs
import com.comsol.model.util.*	% idem
ModelUtil.showProgress(true) ;  % Progress bar
% load all the information needed to re-use the QNM data from previous calculations
load(normalisation_file,'model_name','save_model','sym_factor','omega_var_name',...
    'background_material','total_field_study_name','geometry_name','mesh_name',...
    'electric_point_dipole_name','solver_name','omega','normal_coeff');
pole=omega;
N_coeff=normal_coeff(end);
omegas=2*pi*c./wavelength*1e9;
clear omega
% load the model save after the normalisation process (it contains QNM data)
model = mphload(save_model); 
tic
% Retrieve the domains numbers composing the resonnator (to integrate onto their volumes)
doms=1:model.geom(geometry_name).getNDomains();% list of geometrical domains in the model sheet
doms_bg=model.material(background_material).selection.entities(3);% list of geometrical domains composing the background
doms_res=setdiff(doms,doms_bg);% list of geometrical domains of the resonnator (domains that are not background)

% Retrieve the permitivity (expressions) for all resonator and background materials 
eps_bg_str=char(model.material(background_material).propertyGroup('def').getString('relpermittivity'));
eps_res_str=char(model.material(resonator_material).propertyGroup('def').getString('relpermittivity'));
eps_bg_str=strrep( eps_bg_str,[ '*' omega_var_name '/omega_L_QNM'], '');% to compensate 'setCOMSOL_ComplexFreq'
eps_res_str=strrep( eps_res_str,[ '*' omega_var_name '/omega_L_QNM'], '');% to compensate 'setCOMSOL_ComplexFreq'

% In dataset 1 and 2 (total and background field) we only need data in the volume of the resonator ;
model.result.dataset('dset1').selection.geom(geometry_name,3);
model.result.dataset('dset1').selection.set(doms_res);
model.result.dataset('dset2').selection.geom(geometry_name,3);
model.result.dataset('dset2').selection.set(doms_res);

% Create a join dataset between total field and background field solutions to calculate scattered field
model.result.dataset.create('scatt_field','Join');
model.result.dataset('scatt_field').set('data','dset1');
model.result.dataset('scatt_field').set('data2','dset2');
%by default the join dataset is the difference between data1 and data2

% Create a volume integral operator in COMSOL model sheet 
model.result.numerical.create('int_QNM_PW', 'IntVolume');
model.result.numerical('int_QNM_PW').set('data', 'scatt_field');

% Create a global expression evaluator in model sheet (to evaluate permittivities)
model.result.numerical.create('gev_QNM', 'EvalGlobal');
model.result.numerical('gev_QNM').set('data','dset1');

% update solutions in model sheet (to take  into account the new operators)
model.sol(solver_name).updateSolution();

%% OVERLAP INTEGRAL CALCULATIONS
fprintf('Calculating cross-sections using QNM ...\n');
h=waitbar(0,'Calculate cross-sections with QNM ...');
for i=1:length(wavelength)
    % Retrieve epsilon relative for background
    model.result.numerical('gev_QNM').set('expr',strrep(eps_bg_str,omega_var_name,num2str(omegas(i),'%12.15e')));
    temp=model.result.numerical('gev_QNM').getReal()+1i*model.result.numerical('gev_QNM').getImag();
    eps_bg(i)=mean(temp(:));% materials must be isotropic but it is possible COMSOL gives permittivity tensor (3x3) instead of scalars
    % Retrieve epsilon relative for resonator
    model.result.numerical('gev_QNM').set('expr',strrep(eps_res_str,omega_var_name,num2str(omegas(i),'%12.15e')));
    temp=model.result.numerical('gev_QNM').getReal()+1i*model.result.numerical('gev_QNM').getImag();
    eps_res(i)=mean(temp(:));
    
    % Calculate wave vector for plane wave
    k_diel=sqrt(eps_bg(i))*(omegas(i)./c);% wave vector in dielectric background
    k_diel=num2str(k_diel,'%12.15e');% send it to the model sheet
    
    % Calculate overlap integrals
    model.result.numerical('int_QNM_PW').set('expr',['(emw.Ez)*exp(i*' k_diel '*x)']);%Integral QNM*[PW+k]
    I_overlap_1(i)=N_coeff/sqrt(sym_factor)*(model.result.numerical('int_QNM_PW').getReal()+1i*model.result.numerical('int_QNM_PW').getImag());
    model.result.numerical('int_QNM_PW').set('expr',['(emw.Ez)*exp(-i*' k_diel '*x)']);%Integral QNM*conj([PW+k]) = QNM*[PW-k]
    I_overlap_2(i)=N_coeff/sqrt(sym_factor)*(model.result.numerical('int_QNM_PW').getReal()+1i*model.result.numerical('int_QNM_PW').getImag());
    waitbar(i/length(wavelength),h,'Calculate cross-sections with QNM ...');
end
close(h)
clear temp  h k_diel
if sym_factor >1
    I_overlap=I_overlap_1+I_overlap_2;
    if sym_factor==4
        I_overlap=2*I_overlap;
    end
    I_overlap_conj=I_overlap;% Integral (QNM*PW+k) = Integral (QNM*PW-k) in this particular case with symmetry
else
    I_overlap= I_overlap_1;
    I_overlap_conj= I_overlap_2;
end

% Calculate integral of |QNM|^2
model.result.numerical('int_QNM_PW').set('expr','(emw.normE)^2');
I_QNM2=abs(N_coeff)^2/sym_factor*model.result.numerical('int_QNM_PW').getReal();
I_QNM2=I_QNM2*sym_factor; % |Ez|^2 is integrated here

%% calculate beta and then cross-sections
% coupling coefficient beta of the plane wave to the QNM
beta=-omegas./(omegas-pole).*eps0.*(eps_res-eps_bg).*I_overlap;
% Poynting flux of the plane wave
poynting_PW=1/2*real(sqrt(eps0*eps_bg/mu0)); % Remember; plane wave amplitude is 1 in [V/m]
% integral of |PW|^2
I_PW2=mphint2(model,'1 [V/m]',3,'selection',doms_res)*sym_factor;% plane wave amplitude is 1 in[V/m] 
% calculate the extinction cross-section using Eq.(14) in Qiang Bai's OE paper
sigma_ext=-omegas./poynting_PW/2.*imag( eps0*(eps_res-eps_bg).*(beta.*I_overlap_conj+I_PW2) );
% calculate the absorption cross-section using Eq.(13) in Qiang Bai's OE paper
sigma_abs=-omegas./poynting_PW/2.*imag(eps0*eps_res).*(abs(beta).^2.*I_QNM2+I_PW2+2*real(beta.*I_overlap_conj));% plane wave amplitude is 1[V/m]
% scattering cross section
sigma_scat=sigma_ext-sigma_abs;
toc

h1=figure;
subplot(1,2,1)
plot(wavelength,sigma_abs,'m--');
ylabel('Abs. cross-section (m^2)');
xlabel('wavelength (nm)');
subplot(1,2,2)
plot(wavelength,sigma_scat,'m--');
ylabel('Scatt. cross-section (m^2)');
xlabel('wavelength (nm)');

clear i ans model h1 eps_bg_str eps_res_str omegas choice
save(save_file,'normalisation_file','resonator_material','wavelength','sigma_ext','sigma_scat','sigma_abs');


