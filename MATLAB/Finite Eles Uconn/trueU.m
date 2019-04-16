function out=trueU(in,n)
% Evaluate true solution at input points "in" for test case "n"

% Author: Jieun Lee
% Date: 2/11/2016

if (n==1)
   out= 1 + in;
elseif n==2
   out= 1 + in;
elseif n==3
   out= exp(in) + sin(pi*in);
else
   error('trueU :: invalid test case... exiting.');
end

end % end of function 
