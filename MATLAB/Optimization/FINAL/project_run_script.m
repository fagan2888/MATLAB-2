%Joshua Enxing
%Tufts University
%MA150

%run script for code

m = 10;
n = 2;
eta = 1;

%execute code written by Professor Hu
[A,b] = generate_Ab(m,n,eta);

%build cell array of Hessian functions
hess_f = cell(m,1);
for i=1:m
    hess_f{i} = @(x)0.1*A{i};
end

%build cell array of gradient functions
grad_f = grad_eval(A,b,m);

%parameters for the methods
R = 1;
max_iters = 60;
q = 1/3;
tol = 1E-6;
x0 = [0 0]';
s = 3/4;
mu_0 = 0;
H_0 = 0;
seed = 20;

%use these to find the true minimizer
M = zeros(2);
v = zeros(2,1);

for i=1:10
    M = M + A{i};
    v = v + b{i};
end

%true minimizer
sol = M\v;

%make matrix for k*errors of methods
M = zeros(2,10);
z=max_iters;

%iterate for different values of max_iters
for i = 1:10
    max_iters = z*i;
    [x_sgd,err1] = sgd(grad_f,m,x0,s,max_iters,sol,tol,seed,R);
    M(1,i) = max_iters*norm(x_sgd-sol);
    
    [x_bias,err2,iters2] = bias_removal(grad_f,hess_f,m,x0,s,q,max_iters,sol,tol,mu_0,H_0,seed,R);
    M(2,i) = max_iters*norm(x_bias-sol);
end

%get vector of iteration numbers so can plot against error
iter = zeros(max_iters,1);
for i=1:max_iters
    iter(i) = i;
end

%plot abs error for standard SGD
figure;
plot(iter,err1,iter,1./iter,iter,1./(iter.^s));
legend('abs error','1/k','1/k^s');
title('Absolute Error of Standard SGD');
xlabel('Iteration');
ylabel('||x_k - x^*||');

%plot abs error for BIRR (as called in the paper)
figure;
plot(iters2,err2,iters2,1./iters2);
legend('abs error','1/k');
title('Absolute Error of Bias Removal with q-suffix Averaging');
xlabel('Iteration');
ylabel('||x_k - x^*||');
