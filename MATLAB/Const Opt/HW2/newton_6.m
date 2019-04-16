%Newton's Method for Problem 6, x_0 column vector
function x = newton_6(x_0,tol)

y = x_0;

%Function for x
f = @(x) [x(1).^2+x(2).^2-1 5*x(1).^2-x(2)-2]';

%Function for the Jacobian of f
Jf = @(x) [2*x(1) 2*x(2); 10*x(1) -1];

i = 1;
%Newton iterations
while norm(f(y)) > tol
    y = y - Jf(y)\f(y);
    i = i+1;
end

x = y;