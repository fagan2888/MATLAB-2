function F=Load_New_Old(T,fnf,fnk,g,h,t)

if nargin<5
   h=[];
end
if nargin<4
   g=[];
end
if nargin<3
   fnk=[];
end

% Allocate the load vector F:

Nf=length(T.FNodePtrs);
F=zeros(Nf,1);

% Add the contribution from each element

vo=ones(3,1);
Nt=size(T.Elements,1);
for k=1:Nt

   % Get the coordinates of the vertices and their indices
   % in FNodePtrs and CNodePtrs:

   [c,ll]=getNodes1(T,k);

   if ~isempty(fnf)

      % The lone quadrature point is the centroid of the triangle:

      qpt=(1/3)*sum(c);

      % Compute the area of the triangle:

      J=[c(2,1)-c(1,1) c(3,1)-c(1,1);c(2,2)-c(1,2) c(3,2)-c(1,2)];
      A=0.5*abs(det(J));

      % Get the value of fnf at the quadrature node:

      fval=feval(fnf,qpt(1),qpt(2),t);

      % The triangle contributes to at most 3 entries in the load vector.
      % All contributions are the same (all three basis functions have
      % value 1/3 at the quadrature node):

      I=(A*fval)/3;

      for r=1:3

         % If this vertex is free, then there is a contribution to the
         % load vector:

         llr=ll(r);
         if llr>0

            F(llr)=F(llr)+I;

         end

      end

   end

   % Add the contribution from the Dirichlet data to the right hand side
   % (if the Dirichlet condition is not homogeneous).

   if ~isempty(g)

      % We have nothing to do unless at least one of the vertices of
      % this triangle is a constrained node.

      if any(ll<0)

         % The quadrature point and triangle area must be computed
         % if they were not done above:

         if isempty(fnf)

            % The lone quadrature point is the centroid of the triangle:

            qpt=(1/3)*sum(c);

            % Compute the area of the triangle:

            J=[c(2,1)-c(1,1) c(3,1)-c(1,1);c(2,2)-c(1,2) c(3,2)-c(1,2)];
            A=0.5*abs(det(J));

         end

         % Get the nodal values of the piecewise linear function G
         % agreeing with the Dirichlet data at the constrained nodes
         % and having value zero at each free node.

         w=zeros(3,1);
         for j=1:3
            if ll(j)<0
               w(j)=g(-ll(j));
            end
         end

         % The inverse of M contains the coefficients of the basis
         % functions on this triangle:

         M=[vo,c];
         C=inv(M);

         % Compute the dot products of the gradient of G with the
         % the gradients of the three basis functions.

         h1=C(2:3,:)'*(C(2:3,:)*w);

         % Compute all possible integrals over Tk.  We use a simple
         % one-point formula:

         if isnumeric(fnk)
            I=(A*fnk)*h1;
         elseif ~isempty(fnk)
            I=(A*feval(fnk,qpt(1),qpt(2)))*h1;
         else
            I=A*h1;
         end

         G = zeros(3,3);
         for i=1:3
            for j=1:3
                G(j,i) = ((1/2)*C(2,i)+(1/2)*C(3,i))*(C(1,j)+qpt(1)*C(2,j)+qpt(2)*C(3,j));
            end
         end
         
         % Now add the contributions to the right hand side

         for s=1:3
            lls = ll(s);
            for r = 1:3
                llr=ll(r);
                if llr > 0
                    F(llr)=F(llr)-I(r);
                    if lls < 0
                        F(llr) = F(llr)-G(s,r)*I(r);
                    end
                end
                
            end

         end

      end

   end

   % Now handle the inhomogeneous Neumann data, if any

   if ~isempty(h)

      % Loop over the three edges and compute the contribution
      % of each:

      for j=1:3

         % There's nothing to do unless this edge is free:

         edge=abs(T.Elements(k,j));
         eptr=-T.EdgeEls(edge,2);
         if eptr>0

            % Extract the indices of the endpoints in Nodes:

            ii=T.Edges(edge,1:2);

            % Get the coordinates of the endpoints of the edge:

            c=T.Nodes(ii,1:2);

            % Get the indices of the endpoints in the list of
            % free nodes:

            ll=T.NodePtrs(ii);

            % Evaluate h at the midpoint (use the linear interpolant):

            hval=0.5*sum(h(eptr,:));

            % Evaluate the length of the edge:

            len=norm(c(1,:)-c(2,:));

            % Evaluate the integral (use the midpoint rule):

            I=0.5*len*hval;

            % Now add the contributions to F for each free endpoint:

            for r=1:2
               llr=ll(r);
               if llr>0
                  F(llr)=F(llr)+I;
               end
            end

         end

      end

   end

end
