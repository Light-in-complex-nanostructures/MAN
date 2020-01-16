function [ Y ] = hankel1sph( X, n )
% Spherical Hankel function of the first kind

    Y = sqrt(pi./(2*X)).*besselh(n+0.5, 1, X);

end
