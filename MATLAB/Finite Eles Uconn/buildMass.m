function M=buildMass(grid,dof,element,nElements,nFree)
% Populate the standard mass matrix for 1D CG FEM
% INPUT grid      : mesh nodes 
% INPUT dof       : map global node index to free/constrained node index
% INPUT element   : nodes for each element
% INPUT nElements : number of mesh elements 
% INPUT nFree     : number of free nodes 
% OUTPUT M        : mass matrix 

% Author: Jieun Lee
% Date: 2/09/2016

% initialize output 
M = sparse(nFree, nFree);  

% quadrature nodes
quadNodes= 2;  

for i=1:nElements

  [x, w]= gl_weight(grid(i), grid(i+1), quadNodes);  % Get quadrature weights and abscissa 
  nodes= [grid(i), grid(i+1)];  % Array of nodes for current element (endpoints for first homework) 

    for j=1:2    % Loop over local basis functions 
        for k=1:2

            n_ij= dof(i+j-1);  % Index of local basis function "j" on element "i" within 
            % free or constrained data arrays 
            n_ik= dof(i+k-1);  % Similarly for local basis function "k" on ... 
            if ((n_ij>0)&&(n_ik>0)) 
		    % this excludes the Dirichlet boundary nodes 
            phi_j=evalBasis(x, quadNodes, nodes, 2, j);  
            phi_k=evalBasis(x, quadNodes, nodes, 2, k);
            M(n_ij, n_ik)=M(n_ij, n_ik)+sum(w.*phi_j.*phi_k); 
            end

        end % k-loop over local nodes 
    end % j-loop over local nodes 
end % loop over elements 

end % end function
