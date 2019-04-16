function out=evalBasis(points,numPoints,nodes,numNodes,n)
% evaluate nodal basis function 
% INPUT points    : array of evaluation points -> quadrature abscissa -> x
% from gl_weight
% INPUT numPoints : number of evaluation points 
% INPUT nodes     : array of node locations on an element 
% INPUT numNodes  : number of nodes
% INPUT n         : number of local basis function (left=1, right=2)

% Author: Jieun Lee
% Date: 2/10/2016

if (n>numNodes)  % It is good to try to anticipate ways a piece of code might 
	         % be incorrectly used elsewhere and give yourself a message to                  % help find the bug
error('evalBasis :: invalid basis index... exiting.');
end
if (n<1)
error('evalBasis :: invalid basis index... exiting.');
end

out=ones(numPoints,1);

for i=1:numPoints
	if n==1 % left basis function 
           out(i) = (points(i)-nodes(n+1))/(nodes(n)-nodes(n+1));
    else    
           out(i) = (points(i)-nodes(n-1))/(nodes(n)-nodes(n-1));
	end
end

end % end of function 
