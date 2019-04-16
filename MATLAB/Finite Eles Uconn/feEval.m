function out=feEval(x,grid,dof,element,n,u)
% Evaluate finite element function 
% INPUT x        : array of points on element number `n' where the FE 
%                  function will be evaluated 
% INPUT grid     : node positions 
% INPUT dof      : map from global node index to free/constrained node index
% INPUT element  : element nodes 
% INPUT n        : index of the element where the function will be evaluated
% INPUT u        : FE function nodal values 

% Author: Jieun Lee
% Date: 

nn=length(x);
out=0*x;

% check validity of input data
%if (% LEFT-MOST NODE>min(x))
%error('feEval :: domain of evaluation must coincide with selected element index... exiting.');
%end
%if (% RIGHT-MOST NODE<max(x))
%error('feEval :: domain of evaluation must coincide with selected element index... exiting.');
%end

for j=1:2
    dum=evalBasis(x, 2, [grid(n), grid(n+1)], 2, j);
    out=out+u(n+(j-1))*dum;
end  % loop over basis functions

end % end of function
