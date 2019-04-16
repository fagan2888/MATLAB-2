%Joshua Enxing
%Tufts University
%MA150

%Quadratic least squares with alpha_k chosen by the armijo 
%Stop when abs(0.5*norm(res_current)^2-0.5*norm(res_previous)) < tol
function [a,b,c] = quadratic_least_squares(a_0,b_0,c_0,alpha_0,rho,cc,max_iters,tol)

alpha = alpha_0;
a = a_0;
b = b_0;
c = c_0;
[x,y] = generate_points;
[r,J] = residual_and_jacobian(x,y);
M = (J')*J;

iters = 0;
a_last = 0;
b_last = 0;
c_last = 0;
disp("Iteration || 0.5*norm(r)^2");
for i=1:max_iters
    disp(iters + " || " + 0.5*norm(r(a,b,c))^2);
    if abs(0.5*norm(r(a,b,c))^2-0.5*norm(r(a_last,b_last,c_last))^2) < tol
        break;
    end
    
    pk = -M\((J')*r(a,b,c));
    
    while 0.5*norm(r(a+alpha*pk(1),b+alpha*pk(2),c+alpha*pk(3)))^2 > 0.5*norm(r(a,b,c))^2 + cc*alpha*(((J')*r(a,b,c))')*pk
        alpha = rho*alpha;
    end
    
    a_last = a;
    b_last = b;
    c_last = c;
    
    a = a + alpha*pk(1);
    b = b + alpha*pk(2);
    c = c + alpha*pk(3);
    iters = iters + 1;
end
