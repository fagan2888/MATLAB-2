function K=buildStiffness(grid,dof,element,nElements,nFree,testCase)
% Populate the weighted stiffness matrix for 1D CG FEM
% INPUT grid      : mesh nodes
% INPUT dof       : map global indices to free/constrained indices 
% INPUT element   : nodes for each element 
% INPUT nElements : number of finite elements 
% INPUT nFree     : number of free nodes 
% INPUT testCase  : which test case to use (used in a.m)
% OUTPUT K        : stiffness matrix 

% Author: Jieun Lee
% Date: 2/10/2016

% initialize output 
K = sparse(nFree, nFree); 

% order of quadrature 
quadOrder= 2; 

for i=1:nElements  % Loop over elements

  [x, w]=gl_weight(grid(i), grid(i+1), quadOrder);
  nodes=[grid(i), grid(i+1)];  % Mesh nodes for current mesh element 
  aVals=a(x,testCase);

    for j=1:2
        for k=1:2 

            n_ij= dof(i+j-1);  % Map local node "j" on element "i" to index in 
	        % array of free or constrained data 
            n_ik= dof(i+k-1);  % Map local node "k" on element "i" to... 

            if ((n_ij>0)&&(n_ik>0)) 
		    % this excludes the Dirichlet boundary node 
            phi_j=evalBasisDerivative(x, quadOrder, nodes, 2, j);
            phi_k=evalBasisDerivative(x, quadOrder, nodes, 2, k);
            K(n_ij, n_ik)=K(n_ij, n_ik)+sum(aVals.*w.*phi_j.*phi_k);
            end

        end % k-loop over local nodes 
    end % j-loop over local nodes 
end % loop over elements 

end % end function
