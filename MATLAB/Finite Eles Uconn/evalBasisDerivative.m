function out=evalBasisDerivative(points,numPoints,nodes,numNodes,n)
% evaluate nodal basis function derivative 
% INPUT points    : array of evaluation points 
% INPUT numPoints : number of evaluation points 
% INPUT nodes     : array of node locations on an element
% INPUT numNodes  : number of nodes
% INPUT n         : number of local basis function (left=1, right=2)

% Author: Jieun Lee 
% Date:  2/11/2016

if (n>numNodes)
error('evalBasisDerivative :: invalid basis index... exiting.');
end
if (n<1)
error('evalBasisDerivative :: invalid basis index... exiting.');
end

out=zeros(numPoints,1);

for i=1:numPoints
    if n==1
       out(i) = -1/(nodes(n+1)-nodes(n));
    else
       out(i) = 1/(nodes(n)-nodes(n-1));
    end
end

end % end of function 
