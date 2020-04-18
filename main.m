function [v, err,count] = main (P,r,T,a,b)
% The function is to find roots using the
% newton-raphson method.
% Input P = input pressure.
% Input R is the R constant
% Input T is the temprature value
% Input a and b are the constant values
% Output V = V calculated using the rKwong method
% Output err = magnitude of function
% output count = number of iterations requried
% -------------------------------------------------
% version 1.0: created 15/02/2020 Author: M.Syed
% this matlab function is not flexible, only works
% for one specific function
%---------------------------------------------------
if (P <= 0) | (~isscalar(P)) | (~isreal(P))
error('value P must be a real positive scalar'); % error checking code
end
%--------------------Defining Constants------------------------------%
itt_limit = 3000;
tolerance = 10^-7;
%---------------------Writing the formula------------------------------%
A = a.*P/((r^2)*T^(5/2));
B = (b.*P)/(r*T);
z = 1; %inital estimate
for count=1:itt_limit
if count == itt_limit
error('Iteration limit reached, Iteration did not converge');
break;
end
f = polyval([1 , -1 , (A-B-(B.^2) , -A*B],z) %f(z)
if abs(f) < tolerance
break;
end
df = ([3,-2, (A-B-(B.^2),z] %differential of f(z)
z = z – f./df; %newton-raphson algorithim
v = (z*r*T)./P; %working out the volume.
end
err = abs(f);