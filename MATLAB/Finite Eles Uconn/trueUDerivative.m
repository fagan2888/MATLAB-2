function out=trueUDerivative(in,n)
% Evaluate true solution derivative at input points "in" for test case "n"

% Author: Jieun Lee
% Date: 2/11/2016

if (n==1)
   out= 1;
elseif n==2
   out= 1;
elseif n==3
   out= exp(in) + pi*cos(pi*in);
else
   error('trueU :: invalid test case... exiting.');
end

end % end of function 
