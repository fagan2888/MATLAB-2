%Runs two_grid_tridiag.m max_iters times
function [x,error,res,rel_res] = run_two_grid_tridiag(x_h,b_h,A_h,x_tru,nu_1,max_iters)

x = x_h;

error = zeros(max_iters,1);
res = zeros(max_iters,1);
rel_res = zeros(max_iters-1,1);

for i=1:max_iters
    x = two_grid_tridiag(A_h,x,b_h,nu_1);
    error(i) = norm(x - x_tru);
    res(i) = norm(b_h - A_h*x);
    
    if i > 1
        rel_res(i) = res(i)/res(i-1);
    end
end