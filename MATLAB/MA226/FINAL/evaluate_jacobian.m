function jac = evaluate_jacobian(p)

mat = load('vals.mat');
t = mat.t;

jac = zeros(length(t),4);

jac(:,1) = exp(p(3).*t);
jac(:,2) = exp(p(4).*t);
jac(:,3) = p(1).*t.*exp(p(3).*t);
jac(:,4) = p(2).*t.*exp(p(4).*t);

