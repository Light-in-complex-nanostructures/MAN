function [ y ] = kroneckerDelta_wt( x1, x2 )

if(x1==x2)
    y=1;
else
    y=0;
end

end

