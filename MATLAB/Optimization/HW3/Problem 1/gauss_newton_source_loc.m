%Joshua Enxing
%Tufts University
%MA150

%Nonlinear least squares with alpha_k chosen by the armijo 
%Stop when abs(0.5*norm(res_current)^2-0.5*norm(res_previous)) < tol
function [rs,x] = gauss_newton_source_loc(x_0,d,A,alpha_0,rho,cc,max_iters,tol)

x = x_0;
rx = residual_source_loc(x,d,A);
J = jac_source_loc(x,A);
rs = zeros(max_iters,1);

iters = 0;
disp("Iteration || norm(J'*rx) || f(x)");
for i=1:max_iters
    disp(iters + " || " + norm(J*rx)^2 + " || " + norm(rx)^2);
    
    J = jac_source_loc(x,A);
    if i > 1
       if norm(J*rx)^2 < tol
           rs = rs(1:i-1);
           break;
       end
    end
    
    M = J*(J');
    pk = -M\(J*rx);
    
    alpha = alpha_0;
    rpk = residual_source_loc(x+alpha*pk,d,A);
    while 0.5*norm(rpk)^2 > 0.5*norm(rx)^2 + cc*alpha*((J*rx)')*pk
        alpha = rho*alpha;
        rpk = residual_source_loc(x+alpha*pk,d,A);
    end
    
    x = x + alpha*pk;
    rx = residual_source_loc(x,d,A);
    rs(i) = norm(rx);
    
    iters = iters + 1;
end