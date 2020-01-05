%Steepest descent with sufficient decrease condition for the function
%jhat from HW3. 

%Inputs: initial guess z_0, beta as in suff dec cond, alpha as in suff dec
%cond, max allowable iterations max_iters, and tolerance tol

%Suggested inputs: z_0 = [1 1]', beta = 0.9, alpha = 1E-4, max_iters = 1000,
%tol = 1E-6

%Outputs: computed solution z (which is the desired khat),
%and total number of iterations iters
function [z,iters] = steep_dc_jhat(z_0, beta, alpha, max_iters, tol)

iters = 0;
z = z_0;

%calculation of jhat from hw3
f = @(k) jhat(k);

%calculation of gradient of jhat from hw3
df = @(k) gradient(k);

while norm(df(z)) > tol && iters < max_iters
    %descent direction
    p_k = -1*df(z);
    
    %sufficient decrease
    fxk = f(z);
    p = 1;
    while f(z+p*p_k) > fxk - alpha*p*(p_k)*(p_k)'
        p = p*beta;
    end
    
    %update guess
    z = z + p*p_k;
    
    iters = iters + 1;
end
