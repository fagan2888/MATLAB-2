function A = finite_vol_burger(x,v,T,num_steps)

h_t = T/num_steps;
m = num_steps + 1;
n = length(x) - 2;
h_x = (x(length(x))-x(1))/(n+1);
t = zeros(m,1);

for u=1:m
    t(u) = h_t*(u-1);
end

A = zeros(n+2,m);

A(:,1) = v;
A(1,:) = v(1);
A(n+2,:) = v(n+2);

for j=2:m
    for i = 2:n+1
        A(i,j) = A(i,j-1) - (h_t/(4*h_x))*(A(i+1,j-1)^2 - A(i-1,j-1)^2 - ...
            abs(A(i,j-1)+A(i+1,j-1))*(A(i+1,j-1)-A(i,j-1)) - ...
            abs(A(i,j-1)+A(i-1,j-1))*(A(i-1,j-1)-A(i,j-1)));
    end
end

surf(x(1:n+2),t,transpose(A),'FaceAlpha',0.5);
title('Plot for 1b BCs');
xlabel('x value');
ylabel('time');
zlabel('u value');