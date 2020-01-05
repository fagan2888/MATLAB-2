function [grad,grad_full] = grad_eval(A,b,x,m)

grad = zeros(2,m);
grad_full = zeros(2,1);

for i=1:m
    grad(:,i) = (1/10)*(A{i}*x -b{i});
    grad_full = grad_full + grad(:,i);
end
