function [ rOut, tOut, phiOut ] = Mnm_fixed( n, m, t, phi, r, k )
%Mnm First Spherical Harmonic
% First Spherical Harmonic of degree n and order m calculated at the point
% (r,t,p) in spherical coordinates and given in the spherical base defined
% by t and p (er, et, ep)
% Expression given by S Mühlig & al. Multipole Analysis of meta-atoms

% hankel1sph is the spherical first type hankel bessel function
F = hankel1sph(k*r, n).*exp(1i*m*phi);

rOut = zeros(size(phi));
tOut = F * 1i .* Pinm(cos(t), n, m);
phiOut = -F .* Taunm(cos(t), n, m);

end

