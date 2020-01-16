function [ rOut, tOut, phiOut ] = Nnm_fixed(n, m, t, phi, r, k)
% Nnm Second Spherical Harmonic
% Second Spherical Harmonic of degree n and order m calculated at the point
% (r,t,p) in spherical coordinates and given in the spherical base defined
% by t and p (er, et, ep)
% Expression given by S Muhlig & al. Multipole Analysis of meta-atoms

L = legendre2(n,m,cos(t));
B = hankel1sph(k*r, n) ./ (k*r);
%B = hankel2sph(k*r, n) ./ (k*r);
rOut = n*(n+1) * L .* B .* exp(1i*m*phi);

% derivation 
% h = 1e-9; % at our scale it should be ok to leave this value
% rm = r-h; rp = r+h; dr = 2*h;
% lhv = rp.*sqrt(pi./(2*k*rp)).*besselh(n+0.5, 1, k*rp);
% rhv = rm.*sqrt(pi./(2*k*rm)).*besselh(n+0.5, 1, k*rm);
% derivee = (lhv - rhv)/dr;

% Derivative of r*hankel1sph(k*r,n) with respect to r
 drh = 1/2*sqrt(pi./(2*k*r)).*(besselh(n+0.5, 1, k*r) + k*r.*(besselh(n-0.5, 1, k*r)-besselh(n+1.5, 1, k*r)));
% drh = 1/2*sqrt(pi./(2*k*r)).*(besselh(n+0.5, 2, k*r) + k*r.*(besselh(n-0.5, 2, k*r)-besselh(n+1.5, 2, k*r)));

F = exp(1i*m*phi)./(k*r).*drh;

tOut = Taunm(cos(t),n,m).*F;
phiOut = 1i*Pinm(cos(t),n,m).*F;
end
