function x = interpolate_2D_vec(x)

n = length(x);
m = sqrt(n);
A = vec2mat(x,m);
B = zeros(2*m-1);


for i=1:2*m-1
    for j=1:m
        if mod(i,2)==0
            B(i,2*j-1) = 0.5*(A(i/2,j)+A(i/2+1,j));
        else
            B(i,2*j-1) = A((i+1)/2,j);
        end
    end
end

for k=1:m-1
    B(:,2*k) = 0.5*(B(:,2*k-1)+B(:,2*k+1));
end

x = B(:);

