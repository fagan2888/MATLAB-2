%Joshua Enxing
%Tufts University
%MA150

%Evaluates the symbolic gradient function grad
function g = evaluate_gradient(grad,vars,z)

g = zeros(length(vars),1);
for j=1:length(g)
    g(j) = subs(grad(j),vars,z'); 
end