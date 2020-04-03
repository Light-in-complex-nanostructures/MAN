clear all
clc
tic
addpath(genpath(pwd));

import com.comsol.model.*
import com.comsol.model.util.*

%% COMSOL Livelink setting  

    COMSOL.file='QNMEig_NLnanodisk.mph';                        % Name of COMSOL File
    COMSOL.dataset='dset1';                                                % tag of data set of QNM solution FF
    COMSOL.dataset_SH='dset2';                                             % tag of data set of QNM solution SH
    COMSOL.resonator_domain=[7];                                        % the index of the resonator domain

    model=mphload(COMSOL.file);
    
    Mat.eps_inf=1;                                                                   % Metal permittivity  at omega\to\infty                        
    Mat.eps_b=model.param.evaluate('epsilonb');  % background permittivity  
    Mat.omegap_Res=model.param.evaluate('omegap_Res');% Plasma frequency of metal
    Mat.omega0_Res=model.param.evaluate('omega0_Res');
    Mat.gamma=model.param.evaluate('gamma_Res');                                % damping frequency of meta 
    Mat.eps=@(omega) Mat.eps_inf*(1-Mat.omegap_Res^2./(omega.^2-Mat.omega0_Res.^2-1i*omega*Mat.gamma)); % Metal permittivity function
    Mat.chi2 = 200e-12;                                                                  % Nonlinear tensor element chi 2 (m/V)
    Mat.nsub = model.param.evaluate('nsub');
%-> Physical constants

    Phy.epsilon0=8.854187817e-12;                                                    % vacuum permittivity
    Phy.mu0=4*pi*1e-7;                                                               % vacuum permeability
    Phy.c=1/sqrt(Phy.epsilon0*Phy.mu0);                                              % vacuum light velocity

    Comp.ra = model.param.evaluate('ra');
    Comp.rb = model.param.evaluate('rb');
    Comp.cross_section = pi*Comp.ra*Comp.rb;
    Comp.rs = (sqrt(Mat.eps_b)-Mat.nsub)/(sqrt(Mat.eps_b)+Mat.nsub);

    lambda_target = 1600e-9;
model.param.set('lambda_target',[num2str(lambda_target) '[m]']);
model.geom('geom1').run();
model.mesh('mesh1').run();

model.study('std1').run();
model.study('std2').run();

%% Import modes from COMSOL
%-> QNM angular frequency QNM.omega, QNM normalization factor QNM.QN 
    [QNM.omega,QNM.QN]=mphglobal(model,{'QNM_omega','QN'},'solnum','all','dataset',COMSOL.dataset,'Complexout','on'); 

%-> Read Fields in the sphere sampled at Gaussian point of the mesh elements 

    temp=mpheval(model, 'meshvol', 'solnum', 1,'pattern','gauss','dataset',COMSOL.dataset,'selection', COMSOL.resonator_domain);
    QNM.mesh_vol=temp.d1;   % mesh volume   
    QNM.cord=temp.p;        % mesh coordinates   

    temp=mpheval(model, 'emw.Ex', 'solnum', 'all','pattern','gauss','dataset',COMSOL.dataset,'selection', COMSOL.resonator_domain,'Complexout','on');
    QNM.Ex=temp.d1;         % Electric field x component 

    temp=mpheval(model, 'emw.Ey', 'solnum', 'all','pattern','gauss','dataset',COMSOL.dataset,'selection', COMSOL.resonator_domain,'Complexout','on');
    QNM.Ey=temp.d1;         % Electric field y component 

    temp=mpheval(model, 'emw.Ez', 'solnum', 'all','pattern','gauss','dataset',COMSOL.dataset,'selection', COMSOL.resonator_domain,'Complexout','on');
    QNM.Ez=temp.d1;         % Electric field z component 

%-> QNM angular frequency QNM.omega, QNM normalization factor QNM.QN 
    if strcmp(COMSOL.dataset_SH,COMSOL.dataset); QNM_SH=QNM; else
    [QNM_SH.omega,QNM_SH.QN]=mphglobal(model,{'QNM_omega','QN'},'solnum','all','dataset',COMSOL.dataset_SH,'Complexout','on'); 

    temp=mpheval(model, 'emw.Ex', 'solnum', 'all','pattern','gauss','dataset',COMSOL.dataset_SH,'selection', COMSOL.resonator_domain,'Complexout','on');
    QNM_SH.Ex=temp.d1;         % Electric field x component 

    temp=mpheval(model, 'emw.Ey', 'solnum', 'all','pattern','gauss','dataset',COMSOL.dataset_SH,'selection', COMSOL.resonator_domain,'Complexout','on');
    QNM_SH.Ey=temp.d1;         % Electric field y component 

    temp=mpheval(model, 'emw.Ez', 'solnum', 'all','pattern','gauss','dataset',COMSOL.dataset_SH,'selection', COMSOL.resonator_domain,'Complexout','on');
    QNM_SH.Ez=temp.d1;         % Electric field z component 

    clear temp; 
    end
    QNM.lambda = real(2*pi*Phy.c./QNM.omega)*1e9; %Wavelenght resonant modes [nm]
    QNM.eps=Mat.eps(QNM.omega);

    QNM_SH.lambda = real(2*pi*Phy.c./QNM_SH.omega)*1e9; %Wavelenght resonant modes [nm]
    QNM_SH.eps=Mat.eps(QNM_SH.omega);

    %%% Loading file from COMSOL finished %%%
    
%% -> Computational Setting
Comp.lambda = (1600:1:1750)*1e-9;
Comp.omegav=(2*pi*Phy.c./Comp.lambda);                                      % Frequency range  
Comp.phi = 0;                                                                    % polarization with respect to x axis (phi=0 E along x - phi = pi/2 E along y)
Comp.I0 = 1e13; % Plane wave intensity (W/m2)
Comp.E0 = (sqrt(2*Comp.I0*376.730313667));

NQNM_FF = 4;                                                                % Number of main desired modes at FF for reconstruction
NQNM_SH = 7;                                                                % Number of main desired modes at SH for reconstruction

%% Computing Extinction, absorption cross section

Sol.alpha = zeros(length(QNM.omega),numel(Comp.omegav));
SolSH.alpha = zeros(length(QNM_SH.omega),numel(Comp.omegav));
Sol.ext=zeros(1,numel(Comp.omegav));                                           % Extinction cross section

for iii=1:numel(Comp.omegav)
     
    omega=Comp.omegav(iii);                                                % operating frequency  
    lambda = real(2*pi*Phy.c./omega)*1e9;                                  %Wavelenght [nm]
    
    k=omega/Phy.c*sqrt(Mat.eps_b);                                         % wavenumber in background
    
    p = exp(1i*k*QNM.cord(1,:)*0+1i*k*QNM.cord(2,:)*0+1i*k*QNM.cord(3,:)); % Incident plane wave phase term
    pr = exp(1i*k*QNM.cord(1,:)*0+1i*k*QNM.cord(2,:)*0-1i*k*QNM.cord(3,:)); % Reflected plane wave at z=0 phase term
    
    % Background field definition - valid for a plane wave with k vector
    % along z axis
    E_inc_x=(p+Comp.rs*pr).*Comp.E0*cos(Comp.phi); % Background electric field x-component
    E_inc_y=(p+Comp.rs*pr).*Comp.E0*sin(Comp.phi); % Background electric field y-component  
    E_inc_z=p*0; % Background electric field z-component
    
    S0 = Comp.I0;
   % Field overlap integration between incident field and QNM field 
    E_int=sum( bsxfun(@times,QNM.Ex,E_inc_x.*QNM.mesh_vol)+...
              bsxfun(@times,QNM.Ey,E_inc_y.*QNM.mesh_vol)+...
              bsxfun(@times,QNM.Ez,E_inc_z.*QNM.mesh_vol),...
              2);
  
   % Field overlap integration between complex conjuation of incident field and QNM field    
    E_int_c=sum( bsxfun(@times,QNM.Ex,conj(E_inc_x).*QNM.mesh_vol)+...
                bsxfun(@times,QNM.Ey,conj(E_inc_y).*QNM.mesh_vol)+...
                bsxfun(@times,QNM.Ez,conj(E_inc_z).*QNM.mesh_vol),...
                 2);       
             
   % Prefactor of excitation coefficients 
    
   alpha_QNM=Phy.epsilon0*(-omega*(Mat.eps_inf-Mat.eps_b)-QNM.omega.*(QNM.eps-Mat.eps_inf))...
              ./(omega-QNM.omega);                                          % Excitation coefficient pre-factor

  % used for evaluation of extinction cross section       
    alpha_QNM_current=alpha_QNM*Phy.epsilon0.*(omega*(Mat.eps_inf-Mat.eps_b)-QNM.omega.*(QNM.eps-Mat.eps_inf))...
              ./QNM.QN;                                                     % Extinction cross section pre-factor

  % Extintion cross section
    Sol.ext(iii)=1/2*1/S0*...
                imag(sum(alpha_QNM_current.*E_int_c.*E_int));             % in unit of m^2 (+ sign due to exp(iwt) COMSOL convention)
    Sol.extModes(iii,:)=1/2*1/S0*...
                imag(alpha_QNM_current.*E_int_c.*E_int); 
            
    Sol.alpha(:,iii) = (alpha_QNM.*E_int)./sqrt(QNM.QN); % Modal excitation coefficients
      
    % Current integration for abs cross section
    J_int=abs(Sol.alpha(:,iii)...
        .*-1i*Phy.epsilon0.*(QNM.eps-Mat.eps_inf).*QNM.omega.*...
        sqrt((abs(QNM.Ex./sqrt(QNM.QN)).^2+abs(QNM.Ey./sqrt(QNM.QN)).^2+abs(QNM.Ez./sqrt(QNM.QN)).^2))).^2*QNM.mesh_vol.';
                 
  % Absorption cross section
    Sol.abs(iii)=1/2/S0*...
      (Mat.gamma/(Phy.epsilon0.*Mat.eps_inf.*Mat.omegap_Res.^2)).*sum(J_int); % in unit of m^2
    
    Sol.absModes(iii,:)=1/2/S0*...
      (Mat.gamma/(Phy.epsilon0.*Mat.eps_inf.*Mat.omegap_Res.^2)).*(J_int); % in unit of m^2
  
    Et_x = sum(Sol.alpha(:,iii).*QNM.Ex./sqrt(QNM.QN).*(QNM.eps-Mat.eps_inf)/(Mat.eps(omega)-Mat.eps_inf)); % Total field inside resonator at omega x-component
    Et_y = sum(Sol.alpha(:,iii).*QNM.Ey./sqrt(QNM.QN).*(QNM.eps-Mat.eps_inf)/(Mat.eps(omega)-Mat.eps_inf)); % Total field inside resonator at omega y-component
    Et_z = sum(Sol.alpha(:,iii).*QNM.Ez./sqrt(QNM.QN).*(QNM.eps-Mat.eps_inf)/(Mat.eps(omega)-Mat.eps_inf)); % Total field inside resonator at omega z-component
    
    % Nonlinear polarization vector computation
    % This calculation is valid for a materials with zincblende crystalline
    % symmetry as GaAs/AlGaAs etc...
    % Change this part if a different chi2 tensor is used
    
    PSH_x = 2*Phy.epsilon0*Mat.chi2*(Et_y.*Et_z); % Nonlinear polarization vector at SH x-component
    PSH_y = 2*Phy.epsilon0*Mat.chi2*(Et_z.*Et_x); % Nonlinear polarization vector at SH y-component
    PSH_z = 2*Phy.epsilon0*Mat.chi2*(Et_x.*Et_y); % Nonlinear polarization vector at SH z-component
  
    alphaSH_QNM=(-2*omega)...
               ./(2*omega-QNM_SH.omega);
          
    E_intSH=sum( bsxfun(@times,QNM_SH.Ex,PSH_x.*QNM.mesh_vol)+...
              bsxfun(@times,QNM_SH.Ey,PSH_y.*QNM.mesh_vol)+...
              bsxfun(@times,QNM_SH.Ez,PSH_z.*QNM.mesh_vol),...
              2);
          
    E_intSH_c=sum( bsxfun(@times,QNM_SH.Ex,conj(PSH_x).*QNM.mesh_vol)+...
              bsxfun(@times,QNM_SH.Ey,conj(PSH_y).*QNM.mesh_vol)+...
              bsxfun(@times,QNM_SH.Ez,conj(PSH_z).*QNM.mesh_vol),...
              2);
          
    SolSH.alpha(:,iii) = (alphaSH_QNM.*E_intSH)./sqrt(QNM_SH.QN);
	alphaSH_QNM_current=-2*omega*alphaSH_QNM./QNM_SH.QN;
 
  % Extintion cross section at SH
    SolSH.ext(iii)=1/2*1/S0*...
                imag(sum(alphaSH_QNM_current.*E_intSH.*E_intSH_c));             % in unit of m^2
  % Separate mode contribution to SH extinction  
    SolSH.extModes(iii,:)=1/2*1/S0*...
                imag(alphaSH_QNM_current.*E_intSH.*E_intSH_c); 
         
 % Current integration for abs cross section at SH
     J_intSH=abs(SolSH.alpha(:,iii)...
         .*-1i*Phy.epsilon0.*(QNM_SH.eps-Mat.eps_inf).*QNM_SH.omega.*...
         sqrt((abs(QNM_SH.Ex./sqrt(QNM_SH.QN)).^2+abs(QNM_SH.Ey./sqrt(QNM_SH.QN)).^2+abs(QNM_SH.Ez./sqrt(QNM_SH.QN)).^2))).^2*QNM.mesh_vol.'; %%%
            
%   % Absorption cross section at SH
     SolSH.abs(iii)=1/2/S0*...
       (Mat.gamma/(Phy.epsilon0.*Mat.eps_inf.*Mat.omegap_Res.^2)).*sum(J_intSH); % in unit of m^2
     SolSH.absModes(iii,:)=1/2/S0*...
       (Mat.gamma/(Phy.epsilon0.*Mat.eps_inf.*Mat.omegap_Res.^2)).*(J_intSH); % in unit of m^2

end

zeta = compute_zeta(QNM,QNM_SH,Phy,Mat);

[maxFF, Sol.imax] = maxk(max(abs(Sol.alpha).^2,[],2),NQNM_FF);
[maxSH, SolSH.imax] = maxk(max(abs(SolSH.alpha).^2,[],2),NQNM_SH);

colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840; 0 0.6 0.3; 1 0.2 0.2];

%% Plot cross sections

% Linear extinction cross section

figure('Name','Linear Extinction')
hold on
box on
a = area(2*pi*Phy.c./Comp.omegav*1e9, (Sol.ext)/Comp.cross_section, 'FaceColor', [0 0.25 0.25],'FaceAlpha' , 0.1);
xlabel('Wavelength \lambda_{FF} [nm]'),ylabel('Extinction \sigma_{ext}/\sigma_{geom}')
xlim([min(Comp.lambda)*1e9 max(Comp.lambda)*1e9])
L{1} = ['Reconstruction',num2str(numel(QNM.omega)) ,'QNMs'];

for i=1:length(Sol.imax)
    plot(Comp.lambda*1e9, (Sol.extModes(:,Sol.imax(i)))/Comp.cross_section,'--','linewidth', 1,'Color',colors(i,:));
    L{i+1} = ['Mode ' num2str(i) '- \lambda = ' num2str(floor(QNM.lambda(Sol.imax(i)))) 'nm'];
end
legend(L)
set(gca, 'FontName', 'Times', 'FontSize', 12)
hold off
clear L

% SH cross sections
figure('Name','SH extinction')
hold on
box on
plot(Comp.lambda*1e9, SolSH.ext/Comp.cross_section, 'k', 'linewidth', 1); xlabel('Wavelenght FF \lambda [nm]'),ylabel('\sigma^{(2)}_{ext}/\sigma_{geom}')
L{1} = ['Reconstruction ',num2str(numel(QNM_SH.omega)),' QNMs'];
for i=1:length(SolSH.imax)
    plot(Comp.lambda*1e9, SolSH.extModes(:,SolSH.imax(i))/Comp.cross_section,'--','linewidth', 1,'Color',colors(i,:));
    L{i+1} = ['Mode ' num2str(i) '- \lambda = ' num2str(floor(QNM_SH.lambda(SolSH.imax(i)))) 'nm'];
end
set(gca, 'FontName', 'Times', 'FontSize', 12)
ylim([-max(SolSH.ext)/Comp.cross_section*0.1 max(SolSH.ext)/Comp.cross_section*1.1])
legend(L)
clear L
hold off

figure('Name','SH absorption')

subplot(2,1,1)
hold on
box on
set(gca, 'FontName', 'Times', 'FontSize', 12)
plot(Comp.lambda*1e9, abs(SolSH.abs)/Comp.cross_section, 'k', 'linewidth', 1); xlabel('Wavelenght FF \lambda [nm]'),ylabel('\sigma^{(2)}_{abs}/\sigma_{geom}')
L{1} = 'SH Absorption';
ylim([-max(SolSH.ext)/Comp.cross_section*0.1 max(SolSH.ext)/Comp.cross_section*1.1])
title('SH absoprtion cross section')
legend(L)
clear L
hold off

subplot(2,1,2)
hold on
box on
set(gca, 'FontName', 'Times', 'FontSize', 12)
plot(Comp.lambda*1e9, (SolSH.ext-SolSH.abs)/Comp.cross_section, 'k', 'linewidth', 1); xlabel('Wavelenght FF \lambda [nm]'),ylabel('\sigma^{(2)}_{sca}/\sigma_{geom}')
L{1} = ['Reconstruction ',num2str(numel(QNM_SH.omega)),' QNMs'];
for i=1:length(SolSH.imax)
    plot(Comp.lambda*1e9, (SolSH.extModes(:,SolSH.imax(i))-SolSH.absModes(:,SolSH.imax(i)))/Comp.cross_section,'--','linewidth', 1,'Color',colors(i,:));
    L{i+1} = ['Mode ' num2str(i) '- \lambda = ' num2str(floor(QNM_SH.lambda(SolSH.imax(i)))) 'nm'];
end

ylim([-max(SolSH.ext)/Comp.cross_section*0.1 max(SolSH.ext)/Comp.cross_section*1.1])
title('SH scattering cross section')
legend(L)
clear L
hold off

%% Excitation coefficients

plot_coefficients(Comp,Sol,SolSH,QNM,QNM_SH,Sol.imax,SolSH.imax)

%% Plot QNMs in complex plane
figure(8)
set(gcf,'Position', [10 10 1200 500])
subplot(1,2,1)
plot_QNMcomplexplane(QNM,mean(abs(Sol.alpha).^2,2)'/max(mean(abs(Sol.alpha).^2,2))*200)
subplot(1,2,2)
plot_QNMcomplexplane(QNM_SH,mean(abs(SolSH.alpha).^2,2)'/max(mean(abs(SolSH.alpha).^2,2))*200)

%% Plot overlap tensor
plot_OverlapTensor(zeta,SolSH.imax)

%% Plot near field at SH
plot_nearfieldSH_reconstructed(Comp.omegav(1),SolSH.imax(1),model,COMSOL,Comp,SolSH,QNM_SH,'XZ',[-300,300,-100,500]*1e-9)

time = toc;
disp(['Execution time: ' num2str(floor(time/3600)) 'h ' num2str(floor(time/60-floor(time/3600)*60)) 'm ' num2str(floor(time-floor(time/60)*60)) 's'])