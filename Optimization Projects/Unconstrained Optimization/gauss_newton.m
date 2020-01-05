%Joshua Enxing
%Tufts University
%MA150

%Nonlinear least squares with alpha_k chosen by the armijo 
%Stop when abs(0.5*norm(res_current)^2-0.5*norm(res_previous)) < tol
function x = gauss_newton_source_loc(x_0,alpha_0,rho,cc,max_iters,tol)

alpha = alpha_0;
x = x_0;
rx = residual_source_loc(x,d,A);

iters = 0;
x_last = 0;
r_last = residual_source_loc(x_last,d,A);
disp("Iteration || 0.5*norm(r)^2");
for i=1:max_iters
    disp(iters + " || " + 0.5*norm(r(a,b,c))^2);
    if i > 1
        if abs(0.5*norm(rx)^2-0.5*norm(r_last)^2) < tol
            break;
        end
    end
    
    pk = -M\((J')*rx);
    
    rpk = residual_source_loc(x+alpha*pk,d,A);
    while norm(rpk)^2 > norm(rx)^2 + cc*alpha*(((J')*rx)')*pk
        alpha = rho*alpha;
        rpk = residual_source_loc(x+alpha*pk,d,A);
    end
    
    r_last = rx;
    x = x + alpha*pk;
    rx = residual_source_loc(x,d,A);
    
    iters = iters + 1;
end