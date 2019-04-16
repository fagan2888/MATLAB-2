%Joshua Enxing
%Modified file of the stiffness matrix for advection term

function A=Advection(T)


Nf=length(T.FNodePtrs);
A=sparse(Nf,Nf);

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

   M=[vo,c];
   C=inv(M);


   % Compute the area of the triangle:

   J=[c(2,1)-c(1,1) c(3,1)-c(1,1);c(2,2)-c(1,2) c(3,2)-c(1,2)];
   Area=0.5*abs(det(J));

   % Apply the quadrature rule:
   
   qpt=(1/3)*sum(c);
   I=Area;
   
   % The typical integral we must compute is
   %
   %      integral over Ti b*(grad phi1)*phi2.
   %
   % The gradients are constant vectors, since the basis functions
   % are linear over this triangle. 
   
   G = zeros(3,3);
   for i=1:3
       for j=1:3
           G(j,i) = ((1/2)*C(2,i)+(1/2)*C(3,i))*(1/3);
       end
   end
   
   % This matrix is no longer symmetric, so we account for that

   for s=1:3

      lls=ll(s);
      if lls>0

         for r=1:3

            % If both vertices are free, then there is a contribution
            % to the advection matrix

            llr=ll(r);
            if llr>0
                  A(llr,lls)=A(llr,lls)+G(r,s)*I;
            end

         end

      end

   end

end

