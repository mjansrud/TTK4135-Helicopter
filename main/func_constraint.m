function [ c, ceq ] = func_constraint( z )
%FUNC_CONSTRAINT Summary of this function goes here
%   Detailed explanation goes here
global N mx

alpha = 0.5; 
beta = 20;
lambda_t = (2*pi)/3;

c_temp = -Inf(N,1);
lambda_i = z(1:mx:N*mx);
ek = z(5:mx:N*mx);

for n = 1:N
    c_temp(n) = alpha*exp(-beta*(lambda_i(n)-lambda_t)^2) - ek(n);
end
ceq = [];
c = c_temp;

end
 