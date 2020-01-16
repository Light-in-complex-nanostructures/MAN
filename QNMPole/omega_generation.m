function [new_omegas_set,new_fields_set]=omega_generation(omegas_set,fields_set,delta)
%% omega_generation
% Create a new pole estimation frequency triplet from field calculations on the previous frequency triplet 
% Details found in the Opt. Exp. paper, appendix 2
%% Input are:
%-omegas_set : set of complex frequency on which we calculated a field. Can be (1x1, 1x2, 1x3 complex) 
%-fields_set : set of complex field calculated on omega_set. Must be same dimension as omega_set
%-delta : relative shift to create a new frequencies for the set when there is not enough points to interpolate (<3)
%% Output are:
%-new_omegas_set : new set of frequency to use for the next pole estimation. The last point is always a new one.
%                  if omega_set was <3 points , the new set contains one  more, else, we changed one of those 3 points
%-new_fields_set :  new set of fields to use for the next pole estimation. The last point is always 0 (not yet calculated).
%                  Same dimension as new_omegas_set

%% Calculations
Z=1./fields_set;
Z=Z*max(fields_set);%to improve numerical conditioning of the linear system
norm_omegas=max(abs(omegas_set));
omegas_set_n=omegas_set/norm_omegas;
if length(omegas_set)==1 && length(fields_set)==1 % we cannot yet estimate a pole
    new_omegas_set=[omegas_set,omegas_set*(1-delta)];
    new_fields_set=[fields_set,0];
    
elseif length(omegas_set)==2 && length(fields_set)==2 % we cannot yet estimate a pole
    new_omegas_set=[omegas_set,omegas_set(1)*(1+delta)];
    new_fields_set=[fields_set,0];
    
elseif length(omegas_set)==3 && length(fields_set)==3 % we can estimate a pole
    %Z(i)=(omega(i)-Pole)/(a*omega(i)+b) (for i=1,2,3) : we can write
    %it as a linear system :
    % z1.w1.a + z1.b + p = w1
    % z2.w2.a + z2.b + p = w2
    % z3.w3.a + z3.b + p = w3
    %And put it into matrices; A*(a b Pole)=B and solve it
    A(1,:)=omegas_set_n.*Z;
    A(2,:)=Z;
    A(3,:)=1;
    B=omegas_set_n;
    X=B/A;
    new_omega=norm_omegas*X(3);
    
    [~,index]=sort(abs(fields_set));
    new_omegas_set=[omegas_set(index(3:-1:2)),new_omega];
    new_fields_set=[fields_set(index(3:-1:2)),0];
else
    fprintf('omega_generation : wrong dimensions for the sets (omega and/or field)');
    return
end

end