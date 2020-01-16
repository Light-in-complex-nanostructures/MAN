clear; clc;

import com.comsol.model.*
import com.comsol.model.util.*

%% Input file name and Computing setting  

COMSOL.file='QNMEig_bowtie.mph';                                       % Name of COMSOL File
COMSOL.dataset='dset1';                                                % tag of data set of QNM solution
COMSOL.resonator_domain=4;                                             % the index of the resonator domain

%-> Material Parameters
Mat.eps_inf=1;                                                                   % Metal permittivity  at omega\to\infty 
Mat.omegap_Ag=1.3659e16/sqrt(Mat.eps_inf);                                       % Plasma frequency of metal
Mat.gamma=0.0023*Mat.omegap_Ag*sqrt(Mat.eps_inf);                                % damping frequency of metal 
Mat.eps=@(omega) Mat.eps_inf*(1-Mat.omegap_Ag^2./(omega.^2-1i*omega*Mat.gamma)); % Metal permittivity function
Mat.eps_b=1;                                                                     % background permittivity 

%-> Physical constants

Phy.epsilon0=8.854187817e-12;                                                    % vacuum permittivity
Phy.mu0=4*pi*1e-7;                                                               % vacuum permeability
Phy.c=1/sqrt(Phy.epsilon0*Phy.mu0);                                              % vacuum light velocity                                              % light velocity in vaccum  


%-> Computational Setting 
Comp.omegav=linspace(0.2,0.65,1000)*1.3659e16;                                  % Frequency range  
Comp.K=[0 0 1];                                                                 % wavevector direction, [1 0 0] represents x dir, ..., therein the norm of this vector should be 1
Comp.E=[1 0 0];                                                                 % electric field vector component
Comp.ext=zeros(1,numel(Comp.omegav));                                           % Extinction cross section

%% Read data from COMSOL file 

model=mphload(COMSOL.file);  

%-> QNM angular frequency QNM.omega, QNM normalization factor QNM.QN 
[QNM.omega,QNM.QN]=mphglobal(model,{'QNM_omega','QN'},'solnum','all','dataset',COMSOL.dataset,'Complexout','on'); 

%-> symetry factor sym_factor
[sym_factor]=mphglobal(model,{'sym_factor'},'solnum','1','dataset',COMSOL.dataset,'Complexout','on'); 

%-> Read Fields in the sphere sampled at Gaussian point of the mesh elements 
temp=mpheval(model, 'meshvol', 'solnum', 1,'pattern','gauss','selection', COMSOL.resonator_domain);
QNM.mesh_vol=temp.d1;   % mesh volume   
QNM.cord=temp.p;        % mesh coordinates   

temp=mpheval(model, 'emw.Ex', 'solnum', 'all','pattern','gauss','selection', COMSOL.resonator_domain,'Complexout','on');
QNM.Ex=temp.d1;         % Electric field x component 

temp=mpheval(model, 'emw.Ey', 'solnum', 'all','pattern','gauss','selection', COMSOL.resonator_domain,'Complexout','on');
QNM.Ey=temp.d1;         % Electric field y component 

temp=mpheval(model, 'emw.Ez', 'solnum', 'all','pattern','gauss','selection', COMSOL.resonator_domain,'Complexout','on');
QNM.Ez=temp.d1;         % Electric field z component 

% temp=mpheval(model, 'emw.Hx', 'solnum', 'all','pattern','gauss','selection', COMSOL.resonator_domain,'Complexout','on');
% QNM.Hx=temp.d1;         % Magnetic field x component 
% 
% temp=mpheval(model, 'emw.Hy', 'solnum', 'all','pattern','gauss','selection', COMSOL.resonator_domain,'Complexout','on');
% QNM.Hy=temp.d1;         % Magnetic field y component 
% 
% temp=mpheval(model, 'emw.Hz', 'solnum', 'all','pattern','gauss','selection', COMSOL.resonator_domain,'Complexout','on');
% QNM.Hz=temp.d1;         % Magnetic field z component 

clear temp; 

QNM.eps=Mat.eps(QNM.omega);      % metal permittivity at QNM eigenfrequencies 



%% Computing Extinction, absorption cross section

for iii=1:numel(Comp.omegav)
     
    omega=Comp.omegav(iii);                                                % operating frequency (in general real number) 
    
    k=omega/Phy.c*sqrt(Mat.eps_b);                                         % wavenumber in background
    
    p=exp(1i*k*QNM.cord(1,:)*Comp.K(1)+...
          1i*k*QNM.cord(2,:)*Comp.K(2)+...
          1i*k*QNM.cord(3,:)*Comp.K(3));                                   % plane wave phase                                 
    
   E_inc_x=p*Comp.E(1); E_inc_y=p*Comp.E(2);  E_inc_z=p*Comp.E(3);         % incident electric fields    
   
 
   
   % Field overlap integration between incident field and QNM field 
   E_int=sym_factor.*sum( bsxfun(@times,QNM.Ex,E_inc_x.*QNM.mesh_vol)+...
              bsxfun(@times,QNM.Ey,E_inc_y.*QNM.mesh_vol)+...
              bsxfun(@times,QNM.Ez,E_inc_z.*QNM.mesh_vol),...
              2);
     
   % Field overlap integration between the complex conjugate of incident field and QNM field 
   E_int_c=sym_factor.*sum( bsxfun(@times,QNM.Ex,conj(E_inc_x).*QNM.mesh_vol)+...
              bsxfun(@times,QNM.Ey,conj(E_inc_y).*QNM.mesh_vol)+...
              bsxfun(@times,QNM.Ez,conj(E_inc_z).*QNM.mesh_vol),...
              2);
   
          
   S0=norm(Comp.E)^2*1/2*sqrt(Mat.eps_b*Phy.epsilon0/Phy.mu0)*omega;                % time-averaged poynting power of the incident plane wave        
             
   % Prefactor of excitation coefficients 
   alpha_QNM=Phy.epsilon0*(-omega*(Mat.eps_inf-Mat.eps_b)-QNM.omega.*(QNM.eps-Mat.eps_inf))...
             ./(omega-QNM.omega)...
             ./QNM.QN;
  
  % used for evaluation of extinction cross section       
  alpha_QNM_current=alpha_QNM*Phy.epsilon0.*(-omega*(Mat.eps_inf-Mat.eps_b)-QNM.omega.*(QNM.eps-Mat.eps_inf));

  % Extintion cross section
  Comp.ext(iii)=1/2*1/S0*omega*...
                imag(sum(alpha_QNM_current.*E_int_c.*E_int));             % in unit of m^2
  
 % Current integration for absorption cross section
 J_int=sym_factor.*abs(alpha_QNM.*E_int...
     .*-1i*Phy.epsilon0.*(QNM.eps-Mat.eps_inf).*QNM.omega.*...
     sqrt((abs(QNM.Ex).^2+abs(QNM.Ey).^2+abs(QNM.Ez).^2))).^2*QNM.mesh_vol.';
       
 % Absorption cross section
 Comp.abs(iii)=1/2/S0*omega*...
    (Mat.gamma/(Phy.epsilon0.*Mat.eps_inf.*Mat.omegap_Ag.^2)).*sum(J_int); % in unit of m^2

end;

%% Plot cross sections

% extinction cross section in reduced units
figure; semilogy(Comp.omegav./Mat.omegap_Ag, Comp.ext/(138^2*1e-18));xlabel('\omega/\omega_p'),ylabel('\sigma_{ext}/\lambda_p^2')

% absorption cross section in reduced units
figure; plot(Comp.omegav./Mat.omegap_Ag, Comp.abs/(138^2*1e-18));xlabel('\omega/\omega_p'),ylabel('\sigma_{abs}/\lambda_p^2')

% scattering cross section in reduced units
figure; plot(Comp.omegav./Mat.omegap_Ag, (Comp.ext-Comp.abs)/(138^2*1e-18));xlabel('\omega/\omega_p'),ylabel('\sigma_{sca}/\lambda_p^2')
