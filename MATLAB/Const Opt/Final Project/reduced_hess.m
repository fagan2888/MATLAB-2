%Finds the reduced hessian for inputs: x current iterate, upper and lower
%bound constraints u and l, and exact hessian of f
function R = reduced_hess(x,eps,u,l,hess_f)

n = length(x);
R = hess_f(x(1),x(2));

for i=1:n
    if (u(i) - x(i) <= eps) || (x(i) - l(i) <= eps)
        R(i,:) = zeros(size(R(i,:)));
        R(:,i) = zeros(size(R(:,i)));
        R(i,i) = 1;
    end
end

