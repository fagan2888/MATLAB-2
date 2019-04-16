function F = make_vector(f,h)

x = linspace(0,1,1/h+1);
y = linspace(1,0,1/h+1);
[X,Y] = meshgrid(x,y);
X = X(:);
Y = Y(:);
DT = delaunayTriangulation(X,Y);

NT = length(DT(:,1));
F = zeros(length(DT.Points),1);

for i=1:NT
    
   %Get indices of nodes 
   ind = DT(i,:);
   
   %Use indices of nodes to get coordinates
   coord = zeros(3,2);
   coord(1,:) = DT.Points(ind(1),:);
   coord(2,:) = DT.Points(ind(2),:);
   coord(3,:) = DT.Points(ind(3),:);
   
   cent = (1/3)*sum(coord);
   
   area = 0.5*h^2;
   
   fval = f(cent(1),cent(2));
   
   %All basis elements have value 1/3 at centroid
   val = (area*fval)/3;
   
   for j=1:3

         % If this vertex is free
         if coord(j,1)*coord(j,2) > 0 && coord(j,1) ~= 1 && coord(j,2) ~= 1

            F(ind(j))=F(ind(j))+val;

         end

   end
   
end

o = (1/h+1);
F(1:o+1) = [];

for m=1:o-3
    F(m*(o-2)+1:m*(o-2)+2) = [];
end
u = length(F(:,1));
F(u-o:u) = [];

% P = DT.Points;
% counter = 0;
% for m=1:length(P)
%     if P(m,1)==0 || P(m,2)==0 || P(m,1)==1 || P(m,2)==1
%         F(m-counter) = [];
%         counter = counter+1;
%     end
% end

   
   