function f = func(y)

rng(317,'v4');
A = randn(2,5); 
x = randn(2,1); 
disp(x);
d = sqrt(sum((A-x*ones(1,5)).^2) )+0.05*randn(1,5);

f = (norm(y-A(:,1))^2-d(1)^2)^2+(norm(y-A(:,2))^2-d(2)^2)^2+(norm(y-A(:,3))^2-d(3)^2)^2+(norm(y-A(:,4))^2-d(4)^2)^2+(norm(y-A(:,5))^2-d(5)^2)^2;