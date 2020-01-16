function [Y] = Pinm(X,n,m)
% Pinm First angle-dependent function of degree n and order m
%   As defined in Chapter 4 of 'Bohren and Huffman - Absorption and 
%   scattering of light by small particles' page 94

assert(abs(m) <= n)

Y = m*legendre2(n,m,X)./sqrt(1-X.^2);

% in case the denominator is null
h=1e-5;
if( m == 1)
        Y( (1-X) < h ) = -1/2*n*(1+n);
        Y( (X+1) < h ) = 1/2*n*(1+n)*(-1)^n;
    elseif( m == -1)
        Y( (1-X) < h ) = -1/2;
        Y( (X+1) < h ) = 1/2*(-1)^n;
    else
        Y( (1-X) < h ) = 0;
        Y( (X+1) < h ) = 0;
    end
end

