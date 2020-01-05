grad_f1 = @(x) x - 1;
grad_f2 = @(x) 2*x + 1;

grad_f = {grad_f1,grad_f2};

hess_f1 = @(x)1;
hess_f2 = @(x)2;

hess_f = {hess_f1,hess_f2};

R = 1;
max_iters = 50;
q = 1/2;
tol = 1E-6;
x0 = 1;
m = 2;
sol = 0;
s = 3/4;
mu_0 = 0;
H_0 = 0;
seed = 20;
c = 1;

iter = zeros(max_iters,1);
for i=1:max_iters
    iter(i) = i;
end

M = zeros(5,10);
z=max_iters;
for i = 1:10
    max_iters = z*4*i;
    [x_sgd,err1] = sgd(grad_f,m,x0,s,max_iters,sol,tol,seed,R);
    M(1,i) = x_sgd;
    
    [x_rr,err2] = rand_reshuffle(grad_f,m,x0,s,max_iters,sol,tol,seed,R);
    M(2,i) = x_rr;
    
    [x_rr_avg,err3] = rand_reshuffle_avg(grad_f,m,x0,s,max_iters,sol,tol,seed,R);
    M(3,i) = x_rr_avg;
    
    [x_rr_qavg,err4,iters] = rand_reshuffle_suffix_avg(grad_f,m,x0,s,q,max_iters,sol,tol,seed,R);
    M(4,i) = x_rr_qavg;
    
    [x_bias,err5,iters2] = bias_removal(grad_f,hess_f,m,x0,s,q,max_iters,sol,tol,mu_0,H_0,seed,R);
    M(5,i) = max_iters*x_bias;
end

% figure;
% plot(iter,err1,iter,c./(iter.^s));
% legend('abs error','c/k^s');
% title('Absolute Error of Standard SGD');
% xlabel('Iteration');
% ylabel('||x_k - x^*||');
% 
% figure;
% plot(iter,err2,iter,c./(iter.^s));
% legend('abs error','c/k^s');
% title('Absolute Error of RR');
% xlabel('Iteration');
% ylabel('||x_k - x^*||');
% 
% figure;
% plot(iter,err3,iter,0.5./(iter.^s));
% legend('abs error','0.5/k^s');
% title('Absolute Error of RR with Averaging');
% xlabel('Iteration');
% ylabel('||x_k - x^*||');
% 
% figure;
% plot(iters,err4,iters,c./(iters.^s));
% legend('abs error','c/k^s');
% title('Absolute Error of RR with q-suffix Averaging');
% xlabel('Iteration');
% ylabel('||x_k - x^*||');
% 
% figure;
% plot(iters2,err5,iters2,c./(iters2.^s));
% legend('abs error','c/k^s');
% title('Absolute Error of Bias Removal with q-suffix Averaging');
% xlabel('Iteration');
% ylabel('||x_k - x^*||');
