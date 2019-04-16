function [r,J] = residual_and_jacobian(x,y)

x = x';
y = y';
r = @(a,b,c) a.*x.^2+b.*x+c - y;

J = zeros(30,3);
J(:,1) = x.^2;
J(:,2) = x;
J(:,3) = ones(30,1);
