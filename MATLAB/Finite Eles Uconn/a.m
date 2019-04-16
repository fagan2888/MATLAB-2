function out=a(in,n)
% simple evaluation routine 
% INPUT in : scalar value (independent variable) 
% INPUT  n : for different cases, below 

% Author: Jieun Lee
% Date: 2/10/2016
% You could also make this an input to a()

if n==1
   out=ones(length(in),1);  
elseif n==2
   out=((in+1).*(2-in)*0.5) ;  
elseif n==3
   out=((in+1).*(2-in)*0.5) ; 
else
   error('a.m :: invalid case value... exiting.');
end

end
