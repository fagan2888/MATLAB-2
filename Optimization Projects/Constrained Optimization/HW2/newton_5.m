%Newton's Method for problem 5
function x = newton_5(x_0,tol,a)

%Option to do second part of the problem without writing new function, if a
%is input as four, then we are doing second part of problem
if nargin < 3
    a = 1;
end

%Make vector to store all iterates
x = zeros(100,1);
x(1) = x_0;

%Function f
f = @(x) (x-2).^4+(x-2).^5;

%Function f'(x)
df = @(x) 4*(x-2).^3 + 5*(x-2).^4;

i = 1;
%Newton iterations
while f(x(i)) > tol
    x(i+1) = x(i) - a*f(x(i))/df(x(i));
    i = i+1;
end

x = x(1:i);

x1 = log(x(1:i-1)-2);
x2 = log(x(2:i)-2);

%Plot log of k+1 error vs log of k error to get approximate rate of
%convergence and rate constant
figure;
plot(x1,x2);
xlabel('log|x_k-x*|');
ylabel('log|x_{k+1}-x*|');
title('Plot of Log Error with Line of Best Fit');