function out=feEvalDerivative(x,grid,dof,element,n,u)
% Evaluate finite element function derivative
% INPUT x        : array of points on element number `n' where the FE 
%                  function derivative will be evaluated 
% INPUT grid     : node positions 
% INPUT dof      : map global node indices to free/constrained node indices
% INPUT element  : structure; contains element information 
% INPUT n        : index of the element where the function will be evaluated
% INPUT u        : FE function nodal values 

% Jeff Connors 
% Dec. 2013

nn=length(x);
out=0*x;

% check validity of input data
%if (% FILL IN>min(x))
%error('feEval :: domain of evaluation must coincide with selected element index... exiting.');
%end
%if (% FILL IN<max(x))
%error('feEval :: domain of evaluation must coincide with selected element index... exiting.');
%end

for j=1:2  % loop over two local basis functions on the element
    dum=evalBasisDerivative(x, 2, [grid(n), grid(n+1)], 2, j);
    out=out+u(n+(j-1))*dum;
end  % loop over basis functions

end % end of function
