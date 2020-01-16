function [ Y ] = Taunm( X, n, m )
%Taunm Second angle-dependent function of degree n and order m
%   As defined in Chapter 4 of 'Bohren and Huffman - Absorption and 
%   scattering of light by small particles' page 94
assert(abs(m) <= n)

% T = acos(X);
% 
% if(length(T) > 1)
%     h = min(diff(T))/100;
% else
%     h = T/100;
% end
% Tm = T-h; Tp = T+h; dT = 2*h;
% 
% Y = (legendre2(n,m, cos(Tp))' - legendre2(n,m, cos(Tm))')./dT;

if n>abs(m)
    Y = ( n*X.*legendre2(n,m,X) - (n+m)*legendre2(n-1,m,X) ) ./ (sqrt(1-X.^2));
else
    Y = ( n*X.*legendre2(n,m,X)) ./ (sqrt(1-X.^2));
end

% in case the denominator is null
h=1e-5;
if( m == 1)
    Y( (1-X) < h ) = -1/2*n*(1+n);
    Y( (X+1) < h ) = -1/2*n*(1+n)*(-1)^n;
elseif( m == -1)
    Y( (1-X) < h ) = 1/2;
    Y( (X+1) < h ) = 1/2*(-1)^n;
else
    Y( (1-X) < h ) = 0;
    Y( (X+1) < h ) = 0;
end    

end

