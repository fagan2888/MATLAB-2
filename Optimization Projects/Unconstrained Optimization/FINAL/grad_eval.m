function grad = grad_eval(A,b,m)

grad = cell(m,1);

for i=1:m
    grad{i} = @(x)(1/10)*(A{i}*x -b{i});
end
