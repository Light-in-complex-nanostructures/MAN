function [ Qxx, Qxy, Qxz, Qyx, Qyy, Qyz, Qzx, Qzy, Qzz] = quadrupole_test( anm_qnm, nmax, k0_qnm, namb )

epsilon0 = 8.85418782e-12;
c=2.99792458e8;           % speed of light
z0=1/(c*epsilon0*namb^2); % 
D0 = 6*sqrt(30*pi)/(c*k0_qnm^2*z0*1i);

fac1=-1i*sqrt(6)/3;  % The matematerial paper is wrong;

faca=(anm_qnm(2,nmax+1+2)+ anm_qnm(2,nmax+1-2))*1i;

Qxx=(faca+fac1*anm_qnm(2,nmax+1))*D0;
Qyy=(-faca+fac1*anm_qnm(2,nmax+1))*D0;
Qzz=(1i*sqrt(6)*anm_qnm(2,nmax+1))*D0*2/3;

Qxy=(anm_qnm(2,nmax+1-2)-anm_qnm(2,nmax+1+2))*D0;
Qyx= Qxy;

patch=-1; % The patch is due to the associated legendre function
Qxz= ((anm_qnm(2,nmax+1-1)-anm_qnm(2,nmax+1+1)))*1i*D0*patch;
Qzx= Qxz;

Qyz= ((anm_qnm(2,nmax+1-1)+anm_qnm(2,nmax+1+1)))*D0*patch;
Qzy= Qyz;


end
%------------------------------------------------
