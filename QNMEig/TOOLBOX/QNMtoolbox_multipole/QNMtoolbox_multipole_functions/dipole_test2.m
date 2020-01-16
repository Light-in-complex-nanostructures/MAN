function [px,py,pz,mx,my,mz]=dipole_test2(anm_qnm,bnm_qnm,nmax,k0_qnm,namb)

%% Calcul du moment dipolaire px, py, pz

fprintf('Calcul des moments dipolaires px, py, pz, mx, my, mz ..');tic

epsilon0 = 8.85418782e-12;
c=2.99792458e8; %speed of light
z0=1/(c*epsilon0*namb^2);% Impédance du vide
C0 = sqrt(6*pi)*1i/(c*k0_qnm*z0);

% note the minus sign is added to make the electic dipole has a radiation
% pattern of 1/4/pi*k0^2*(n X p) X p
fac_patch1=-1;
px=fac_patch1*C0*(anm_qnm(1,nmax+2)- anm_qnm(1,nmax));  
py=fac_patch1*C0*1i*(anm_qnm(1,nmax+2) + anm_qnm(1,nmax));
pz=-C0*sqrt(2)*anm_qnm(1,nmax+1);

% The minus sign in contrast to the Ref. are caused by the different
% defination of Pnm

fac_patch2=1i*namb;
mx=-c*C0*(bnm_qnm(1,nmax+2)- bnm_qnm(1,nmax))/fac_patch2;      % The minus sign is added in 28/08/2018 
my=-c*C0*1i*(bnm_qnm(1,nmax+2) + bnm_qnm(1,nmax))/fac_patch2;  % The minus sign is added in 28/08/2018
mz=-c*C0*sqrt(2)*bnm_qnm(1,nmax+1)/fac_patch2;

fprintf('done !');toc

end
