%Joshua Enxing
%Tufts University
%MA150

%Evaluates the symbolic hessian function hess
function h_mat = evaluate_hessian(hess,vars,z)

h_mat = zeros(length(vars));
for k=1:length(vars)
    for l=1:length(vars)
        h_mat(k,l) = subs(hess(k,l),vars,z');
    end
end