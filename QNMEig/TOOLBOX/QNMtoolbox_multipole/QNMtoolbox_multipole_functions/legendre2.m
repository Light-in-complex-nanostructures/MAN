function [ L ] = legendre2(n,m,X)
% legendre2: Custom Legendre polynom for negative m values 

L = 1;
if(m < 0)
    m = abs(m);
    % factor found not without shame on wikipedia
    L = (-1)^m*factorial(n-m)/factorial(n+m);
    %L = (-1)^m*sqrt(factorial(n-m)/factorial(n+m));
end
L = L*legendre(n,X);

% legendre returns a whole matrix for all legendre polynoms evaluated at X
% for m values going from 0 to n-1, so we take only the one of interest
L = L(m+1, :).';

end

