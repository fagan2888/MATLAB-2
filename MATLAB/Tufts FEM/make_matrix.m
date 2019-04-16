
%make matrix A from input h value
function A = make_matrix(h)

x = linspace(0,1,1/h+1);
y = linspace(1,0,1/h+1);
[X,Y] = meshgrid(x,y);
X = X(:);
Y = Y(:);
DT = delaunayTriangulation(X,Y);

NT = length(DT(:,1));
A = zeros(length(DT.Points),length(DT.Points));

%loop over elements
for i=1:NT
    
   %Get indices of nodes 
   ind = DT(i,:);
   
   %Use indices of nodes to get coordinates
   coord = zeros(3,2);
   coord(1,:) = DT.Points(ind(1),:);
   coord(2,:) = DT.Points(ind(2),:);
   coord(3,:) = DT.Points(ind(3),:);
   
   %Each basis function on this element will be of the form a + bx + cy, 
   %so we find a, b, and c
   M = ones(3);
   M(:,2:3) = coord;
   C = inv(M);
   
   %Since our functions are linear, the gradients will just be (b,c), so
   %the dot products of each of the gradients will be the following matrix
   G = C(2:3,:)'*C(2:3,:);
   
   area = 0.5*h^2;
   
   for j=1:3
       
       %If node is free
       if coord(j,1)*coord(j,2) > 0 && coord(j,1) ~= 1 && coord(j,2) ~= 1 
           
           for k = 1:j
               
               %If node is free
               if coord(k,1)*coord(k,2) > 0 && coord(k,1) ~= 1 && coord(k,2) ~= 1 
   
                    if ind(k) <= ind(j)
                        A(ind(k),ind(j))=A(ind(k),ind(j))+G(k,j)*area;
                    else
                        A(ind(j),ind(k))=A(ind(j),ind(k))+G(k,j)*area;
                    end
                    
               end
               
           end
           
       end
       
   end 
   
end
    
A=A+triu(A,1)';

o = (1/h+1);
A(1:o+1,:) = [];
A(:,1:o+1) = [];

for m=1:o-3
    A(m*(o-2)+1:m*(o-2)+2,:) = [];
    A(:,m*(o-2)+1:m*(o-2)+2) = [];
end
u = length(A(:,1));
A(u-o:u,:) = [];
A(:,u-o:u) = [];
% o = (1/h+1);
% for m=2:o-1
%     A(1+(o-2)*(m-2):(o-2)*(m-1),1+(o-2)*(m-2):(o-2)*(m-1)) = A((m-1)*o+2:m*o-1,(m-1)*o+2:m*o-1);
% end
% 
% A = A(1:(o-2)^2,1:(o-2)^2);
% P = DT.Points;
% counter = 0;
% for m=1:length(P)
%     if P(m,1)==0 || P(m,2)==0 || P(m,1)==1 || P(m,2)==1
%         A(m-counter,:) = [];
%         A(:,m-counter) = [];
%         counter = counter+1;
%     end
% end


