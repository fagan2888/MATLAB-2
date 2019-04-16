%Joshua Enxing
%Code written to approximate the boundary during mesh refinements
%Now handles case when k-1 new points are desired in between p1 and p2

function pts = parabolaf(p1,p2,k)

if nargin<3
   k=2;
end

pts = zeros(k-1,2);

h = abs(p2(1)-p1(1))/k;

for i = 1:k-1
    if (p2(1) > p1(1))
        pts(i,1) = p1(1) + i*h;
        pts(i,2) = (pts(i,1))^2;
        
    else
        pts(i,1) = p2(1) + i*h;
        pts(i,2) = (pts(i,1))^2;
    end
end
