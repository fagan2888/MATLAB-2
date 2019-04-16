%Joshua Enxing
%Modified code given by Gockenbach for stiffness matrix to make
%Mass matrix

function M=Mass(T)

Nf=length(T.FNodePtrs);
M=sparse(Nf,Nf);

% Loop over the elements, adding the contributions from each.

vo=ones(3,1);
Nt=size(T.Elements,1);
for k=1:Nt

   % Get the coordinates of the vertices and their indices
   % in FNodePtrs and CNodePtrs:

   [c,ll]=getNodes1(T,k);

   % Each basis function, restricted to this triangle, is a function
   % of the form z=g(1)+g(2)x+g(3)y.  Compute the vector g for each of
   % three basis functions (the three vectors are stored as the columns
   % of the matrix C).

   R=[vo,c];
   C=inv(R);

   % The typical integral we must compute is
   %
   %      integral over Ti phi1.phi2

   G = zeros(3,3);
   
   %Barycentric Coordinates
   
   v1_x = (1/6)*c(1,1)+(1/6)*c(2,1)+(2/3)*c(3,1);
   v2_x = (1/6)*c(1,1)+(2/3)*c(2,1)+(1/6)*c(3,1);
   v3_x = (2/3)*c(1,1)+(1/6)*c(2,1)+(1/6)*c(3,1);
   
   v1_y = (1/6)*c(1,2)+(1/6)*c(2,2)+(2/3)*c(3,2);
   v2_y = (1/6)*c(1,2)+(2/3)*c(2,2)+(1/6)*c(3,2);
   v3_y = (2/3)*c(1,2)+(1/6)*c(2,2)+(1/6)*c(3,2);
   
   %Compute phi_i*phi_j with weights
   
   for i=1:3
       for j=1:3
           G(i,j) = (1/6)*((C(:,i)'*[1 v1_x v1_y]')*(C(:,j)'*[1 v1_x v1_y]')+... 
           (C(:,i)'*[1 v2_x v2_y]')*(C(:,j)'*[1 v2_x v2_y]') +...
           (C(:,i)'*[1 v3_x v3_y]')*(C(:,j)'*[1 v3_x v3_y]'));
       end
   end

   % Compute the area of the triangle:

   J=[c(2,1)-c(1,1) c(3,1)-c(1,1);c(2,2)-c(1,2) c(3,2)-c(1,2)];
   A=0.5*abs(det(J));

   % The triangle contributes to at most 6 entries in the (upper triangle
   % of the) mass matrix.  We compute these six entries in the
   % following double loop.

   for s=1:3

      lls=ll(s);
      if lls>0

         for r=1:s

            % If both vertices are free, then there is a contribution
            % to the mass matrix

            llr=ll(r);
            if llr>0

               if llr<=lls
                  M(llr,lls)=M(llr,lls)+G(r,s)*A;
               else
                  M(lls,llr)=M(lls,llr)+G(r,s)*A;
               end

            end

         end

      end

   end

end

% Now fill in the lower triangle of Mass, using the symmetry.

M=M+triu(M,1)';
