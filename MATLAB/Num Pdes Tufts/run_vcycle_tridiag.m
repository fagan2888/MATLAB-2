%Runs vcycle_tridiag.m max_iters times
function [x,error,res,rel_res] = run_vcycle_tridiag(x_h,b_h,A_h,x_tru,nu_1,nu_2,max_iters)

x = x_h;
error = zeros(max_iters,1);
res = zeros(max_iters,1);
rel_res = zeros(max_iters-1,1);

for i=1:max_iters
    x = vcycle_tridiag(x,b_h,A_h,nu_1,nu_2);
    error(i) = norm(x-x_tru);
    res(i) = norm(b_h - A_h*x);
    
    if i > 1
        rel_res(i-1) = res(i)/res(i-1);
    end
end