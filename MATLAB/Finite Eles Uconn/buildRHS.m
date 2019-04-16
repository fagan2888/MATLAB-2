function rhs=buildRHS(grid,dof,element,nElements,nFree,bndyData,testCase)
% Populate the right-hand side forcing vector for 1D CG FEM code
% INPUT grid        : node positions
% INPUT dof         : map global node indices to free/constrained indices 
% INPUT element     : nodes for each element 
% INPUT nElements   : number of elements
% INPUT nFree       : number of degrees of freedom 
% INPUT bndyData    : Dirichlet boundary data array 
% INPUT testCase    : problem parameter 
% OUTPUT rhs

% Author: Jieun Lee
% Date: 2/11/2016

rhs=zeros(nFree,1); 
quadOrder=2;

for i=1:nElements


    [x, w]=gl_weight(grid(i), grid(i+1),quadOrder);
    nodes= [grid(i), grid(i+1)];  %Mesh nodes on current mesh element 
    fVals=f(x, testCase);
    aVals=a(x, testCase);

    for j=1:2  % Loop over basis functions on current element 
        
       n_ij= dof(i+j-1);  % Map mesh node for basis function "j" on element "i" 
                 % to corresponding index within free or constrained data array 
	if (n_ij>0)
           % this test excludes the Dirichlet boundary node basis function
	   % from being counted as if for a free node 
	   phi_j=evalBasis(x, quadOrder, nodes, 2, j);
	   rhs(n_ij)=rhs(n_ij)+sum(fVals.*w.*phi_j);
	   for k=1:2  % Search for local basis functions "k" that correspond to 
		      % boundary nodes and add boundary data to load vector 
	   n_ik= dof(i+k-1);  % Map mesh node for basis function "k" on element                               % "i" to corresponding index within free or 
	                      % constrained data array 
	   if (n_ik<0) % test if this node is on the boundary 
              phi_k=evalBasis(x, quadOrder, nodes, 2, k);
	      dphi_j_dx=evalBasisDerivative(x, quadOrder, nodes, 2, j);
	      dphi_k_dx=evalBasisDerivative(x, quadOrder, nodes, 2, k);
              rhs(n_ij)=rhs(n_ij)-...
	      bndyData(-n_ik)*...
	      sum(aVals.*w.*dphi_j_dx.*dphi_k_dx+w.*phi_j.*phi_k);  
	   end % boundary contributions code block 
           end
        end % if test to exclude boundary basis functions as being free

    end  % loop over local nodes

end % loop over elements

end % end of function
