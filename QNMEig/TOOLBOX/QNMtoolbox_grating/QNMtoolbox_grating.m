clear; clc;

import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.setServerBusyHandler(ServerBusyHandler(2));
%% Input file name and Computing setting  
tic
File_COMSOL=char('QNMEig_grating_theta.mph');                    % COMSOL File Name
COMSOL.resonator_domain='sel2';                                    % Comsol Resonator Domain
COMSOL.dataset='dset1';                                            % tag of data set of QNM solution

model=mphload(File_COMSOL);                                        % Load COMSOL File

%% Material parameters 
Mat.omegap=2*pi;                 % reduced plasma angular frequency of metal
Mat.gamma=mphglobal(model,{'gamma'},'solnum','1','dataset',COMSOL.dataset,'Complexout','on');     % reduced damping frequency of metal
Mat.epsiloninf=1;                % metal permittivity at omega-->infinity
Mat.epsilonb=1;                  % Background permittivity
Mat.omega0=mphglobal(model,{'omega0'},'solnum','1','dataset',COMSOL.dataset,'Complexout','on');                    % Lorentz pole angular resonance frequency
Mat.eps=@(omega) Mat.epsiloninf.*(1-Mat.omegap.^2./(omega.^2-Mat.omega0^.2+1i.*Mat.gamma.*omega)); % Metal permittivity 

fp=1.26e16/(2*pi);               % Plasma frequency in Hz
lp_m=3e8/fp;                     % Plasma wavelength in m

%% Geometry parameters 
% Geometry data is extracted directly from COMSOL model

Geo.a=mphglobal(model,{'a'},'solnum','1','dataset',COMSOL.dataset,'Complexout','on'); % grating period in reduced units
Geo.f=mphglobal(model,{'f'},'solnum','1','dataset',COMSOL.dataset,'Complexout','on'); % filling fraction in reduced units
Geo.height=mphglobal(model,{'height'},'solnum','1','dataset',COMSOL.dataset,'Complexout','on'); % groove height in reduced units
Geo.w=Geo.f.*Geo.a;         %length of metal bump
theta=mphglobal(model,{'theta'},'solnum','1','dataset',COMSOL.dataset,'Complexout','on'); % incidence angle
altitude = mphglobal(model,{'altitude'},'solnum','1','dataset',COMSOL.dataset,'Complexout','on');


%% Extract QNM Fields and frequencies from COMSOL model 


% QNM frequencies
[QNM.omega]=mphglobal(model,{'lambda'},'solnum','all','dataset',COMSOL.dataset,'complexout','on');    % read eigenvalue of COMSOL 

% QNM mode Normalization coefficient
[QNM.QN]=mphglobal(model,{'QN'},'solnum','all','dataset',COMSOL.dataset,'complexout','on');    % read eigenvalue of COMSOL 

% Metal permittivity at QNM eigenfrequencies
Mat.eps_QNM=Mat.epsiloninf.*(1-Mat.omegap.^2./(QNM.omega.^2+1i.*Mat.gamma.*QNM.omega));

% H field integral above grating
% [intHzr0]=mphglobal(model,{'intop5(Hz)'},'solnum','all','dataset',COMSOL.dataset,'complexout','on');    % read integral of H over a period in COMSOL 
% intHzr0=intHzr0/Geo.a;

% Extracting mesh coordinates in Gaussian pattern with weights
temp=mpheval(model, 'meshvol', 'solnum', 1,'dataset',COMSOL.dataset,'pattern','gauss','selection',COMSOL.resonator_domain);
QNM.mesh_vol=temp.d1;   % mesh volume   
QNM.cord=temp.p;        % mesh coordinates

%QNM fields in resonator domain in resonator
QNM.Ex=[];
QNM.Ey=[];
QNM.Ex_m=[];
QNM.Ey_m=[];
QNM.Hz=[];

tic

% QNM E field x component
temp=mpheval(model, 'real(Ex)','outersolnum','all', 'solnum', 'all','dataset',COMSOL.dataset,'pattern','gauss','selection', COMSOL.resonator_domain);
temp2=mpheval(model, 'imag(Ex)','outersolnum','all', 'solnum', 'all','dataset',COMSOL.dataset,'pattern','gauss','selection', COMSOL.resonator_domain);

QNM.Ex=[QNM.Ex;temp.d1+1i.*temp2.d1];         % Electric field x component
% QNM E field y component
temp=mpheval(model, 'real(Ey)','outersolnum','all', 'solnum', 'all','dataset',COMSOL.dataset,'pattern','gauss','selection', COMSOL.resonator_domain);
temp2=mpheval(model, 'imag(Ey)','outersolnum','all', 'solnum', 'all','dataset',COMSOL.dataset,'pattern','gauss','selection', COMSOL.resonator_domain);

QNM.Ey=[QNM.Ey;temp.d1+1i.*temp2.d1];         % Electric field y component
%Ex_m - Counterpropagative field x component
temp=mpheval(model, 'real(Ex_m)','outersolnum','all', 'solnum', 'all','dataset',COMSOL.dataset,'pattern','gauss','selection', COMSOL.resonator_domain);
temp2=mpheval(model, 'imag(Ex_m)', 'solnum', 'all','dataset',COMSOL.dataset,'pattern','gauss','selection', COMSOL.resonator_domain);

QNM.Ex_m=[QNM.Ex_m;temp.d1+1i.*temp2.d1];         % Counterpropagating Electric field x component

%Ey_m - Counterpropagative field y component
temp=mpheval(model, 'real(Ey_m)', 'solnum', 'all','dataset',COMSOL.dataset,'pattern','gauss','selection', COMSOL.resonator_domain);
temp2=mpheval(model, 'imag(Ey_m)', 'solnum', 'all','dataset',COMSOL.dataset,'pattern','gauss','selection', COMSOL.resonator_domain);

QNM.Ey_m=[QNM.Ey_m;temp.d1+1i.*temp2.d1];         % Counterpropgating Electric field y component

%QNM H field z component
temp=mpheval(model, 'real(Hz)', 'solnum', 'all','dataset',COMSOL.dataset,'pattern','gauss','selection', COMSOL.resonator_domain);
temp2=mpheval(model, 'imag(Hz)', 'solnum', 'all','dataset',COMSOL.dataset,'pattern','gauss','selection', COMSOL.resonator_domain);

QNM.Hz=[QNM.Hz;temp.d1+1i.*temp2.d1];         % Magnetic field z component


numt=150;
[xss,yss]=meshgrid(linspace(-Geo.a/2,Geo.a/2,numt),linspace(-Geo.height/2-Geo.a/2,1.5,numt));
xsss=reshape(xss,1,numt^2); ysss=reshape(yss,1,numt^2);
[QNM.Hz_QNM]=mphinterp(model,{'Hz'},'solnum','all','coord',[xsss;ysss],'Complexout','on','dataset',COMSOL.dataset);   % read fields on your specified grids
QNM.Hz_QNM=reshape(QNM.Hz_QNM,length(QNM.omega),numt,numt);


% stop 
 %% QNM field at upper plane   
 
 temp_up=mpheval(model, 'meshvol', 'solnum', 1,'dataset',COMSOL.dataset,'pattern','gauss','selection','sel1');
 QNM.mesh_vol_up=temp_up.d1;   % mesh volume
 cord_up=temp_up.p;        % mesh coordinates
 QNM.x_up=cord_up(1,:); % X coordinates in the upper plane
 QNM.y_up=cord_up(2,:); % Y coordiantes in the upper plane
 
 QNM.Ex_up=[];
 QNM.Ey_up=[];
 QNM.Hz_up=[];
  
%  % Electric field x component 
%  temp=mpheval(model, 'Ex', 'solnum', 'all','dataset',COMSOL.dataset,'pattern','gauss','selection', 'sel1','Complexout','on');
%  QNM.Ex_up=[QNM.Ex_up;temp.d1];         % Electric field x component
%  % Electric field x component 
%  temp=mpheval(model, 'Ey', 'solnum', 'all','dataset',COMSOL.dataset,'pattern','gauss','selection', 'sel1','Complexout','on');
%  QNM.Ey_up=[QNM.Ey_up;temp.d1];         % Electric field y component

% Magnetic field z component 
 temp=mpheval(model, 'Hz', 'solnum', 'all','dataset',COMSOL.dataset,'pattern','gauss','selection', 'sel1','Complexout','on');
 QNM.Hz_up=[QNM.Hz_up;temp.d1];         % Electric field y component
 
t_extract_fields=toc
intHz_gauss=QNM.Hz_up./Geo.a*QNM.mesh_vol_up.'; % QNM 0 order Fourier coefficient



%%

lp=149.6;
a = lp./linspace(250,1000,10);a=sort(a); % excitation frequency from UV --> NIR range
a2 = linspace(a(1),a(end),3000);

tic
for k=1:1:length(a2)

omega=Mat.omegap*a2(k); % Working frequency

eta=sqrt(Mat.eps(omega)-Mat.epsilonb*sin(theta)^2)/(sqrt(Mat.epsilonb)*cos(theta))*Mat.epsilonb/Mat.eps(omega);
rfresnel=-(eta-1)/(eta+1);   % reflection, defined with z-component magnetic field
tfresnel=1+rfresnel;         % transmission, defined with z-component magnetic field 
ks=omega*sqrt(Mat.eps(omega)-Mat.epsilonb*sin(theta)^2);   % perpendicular wavevector modulus inside metal 

% This is the field inside the metal slab in the pertubation domain, the phase is 0 at y=0
Ex_inc=tfresnel*ks./(omega*Mat.eps(omega)).*exp(-1i*ks.*(QNM.cord(2,:)));   % Background field with phase corresponding to slab background
Ey_inc=tfresnel.*Mat.epsilonb.*sin(theta)./(Mat.eps(omega)).*exp(-1i.*ks.*(QNM.cord(2,:))); % Background field with phase corresponding to slab background


intEx=QNM.Ex_m.*Ex_inc*QNM.mesh_vol.'; % Overlap integral x component
intEy=QNM.Ey_m.*Ey_inc*QNM.mesh_vol.'; % Overlap integral y component


% Overlap integral in resonator domain above the substrate
intE=intEx+intEy; % Total integral inside resonator domain performed with Gaussian sampling

alpha_QNM=(Mat.eps(omega)-Mat.epsiloninf-QNM.omega.*(Mat.epsilonb-Mat.eps(omega))...
    ./(omega-QNM.omega))./QNM.QN;  % excitation coefficients prefactor

r0=sum(bsxfun(@times,intHz_gauss,alpha_QNM.*intE))+...
              rfresnel*exp(1i*sqrt(Mat.epsilonb)*omega*cos(theta)*(altitude)); % 0 order fourier scattered field + up going bakground field at y = altitude
       
R0(k)=r0*conj(r0); % Specular Reflection coefficient
end
t_compute_spectra=toc
% trace_reflection=figure;
figure;
f1=a2.*fp;
wavelength=3e14./f1;

plot(wavelength,R0,'-r');

xlabel('\lambda(µm)');ylabel('Reflection')

%% Plot scattered Magnetic field field in space around grating

a3=a2(1000) %value of a2 corresponding to an excitation wavelength of 500 nm

omega=Mat.omegap*a3; % Working frequency

eta=sqrt(Mat.eps(omega)-Mat.epsilonb*sin(theta)^2)/(sqrt(Mat.epsilonb)*cos(theta))*Mat.epsilonb/Mat.eps(omega);
rfresnel=-(eta-1)/(eta+1);   % reflection, defined with z-component magnetic field
tfresnel=1+rfresnel;         % transmission, defined with z-component magnetic field 
ks=omega*sqrt(Mat.eps(omega)-Mat.epsilonb*sin(theta)^2);   % perpendicular wavevector modulus inside metal 

% This is the field inside the metal slab in the pertubation domain, the phase is 0 at y=0
Ex_inc=tfresnel*ks./(omega*Mat.eps(omega)).*exp(-1i*ks.*(QNM.cord(2,:)));   % Background field with phase corresponding to slab background
Ey_inc=tfresnel.*Mat.epsilonb.*sin(theta)./(Mat.eps(omega)).*exp(-1i.*ks.*(QNM.cord(2,:))); % Background field with phase corresponding to slab background


intEx=QNM.Ex_m.*Ex_inc*QNM.mesh_vol.'; % Overlap integral x component
intEy=QNM.Ey_m.*Ey_inc*QNM.mesh_vol.'; % Overlap integral y component

% Overlap integral in resonator domain above the substrate
intE=intEx+intEy; % Total integral inside resonator domain performed with Gaussian sampling

alpha_QNM=(Mat.eps(omega)-Mat.epsiloninf-QNM.omega.*(Mat.epsilonb-Mat.eps(omega))...
    ./(omega-QNM.omega))./QNM.QN; 
    
Hz_reconst=sum(alpha_QNM.*intE.*QNM.Hz_QNM,1); % reconstructed scattered field
Hz_reconst=reshape(Hz_reconst,numt,numt); 

figure,imagesc(xsss,ysss,abs(Hz_reconst.*exp(i.*omega.*(xss-Geo.a/2))));
set(gca,'YDir','normal') 
colorbar;
colormap('hot')
axis equal;
xlabel('x');
ylabel('y');

hold on
% Draw in the edges of the metal 
x=[-Geo.a/2, -Geo.a/2*(0.9) -Geo.a/2*(0.9) +Geo.a/2*(0.9) +Geo.a/2*(0.9) Geo.a/2];
y=[-Geo.height -Geo.height 0 0 -Geo.height -Geo.height];
plot(x,y,'b','linewidth',2);


