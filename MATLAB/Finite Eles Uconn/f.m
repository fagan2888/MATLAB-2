function out=f(in,n)
% simple evaluation rountine for RHS forcing of PDE
% INPUT in : scalar value (independent variable) 
% INPUT  n : for different cases, below 

% Author: Jieun Lee
% Date: 2/11/2016
if n==1
   out= 1 + in;
elseif n==2
   out= 2*in + 0.5;
elseif n==3
   out=0.5.*((2.*in-1).*(exp(in)+pi.*cos(pi.*in))+(in.^2-in-2).*(exp(in)-(pi.^2).*sin(pi.*in)))+exp(in)+sin(pi.*in);
else
   error('f.m :: invalid case value... exiting.');
end

end
