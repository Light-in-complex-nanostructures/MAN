function [epsln] = Drude_epsln_multi( omega, material, imat)
% The subroutine calculate the permitivity
% This subroutine assume that the time dependent part is e[-iwt]
% Which is in contrast with comsol

eps_inf=conj(material(imat).epsln_inf);  % The conjugate is due to Comsol
gamma=material(imat).gamma;
omega0=material(imat).omega0;
omega_p=material(imat).omega_p;
npole=material(imat).npole;

epsln=eps_inf;
for ipole=1: npole
    dno=omega.^2-omega0(ipole).^2+1i.*omega.*gamma(ipole);
    epsln=epsln-eps_inf*omega_p(ipole).^2./dno;
end

end