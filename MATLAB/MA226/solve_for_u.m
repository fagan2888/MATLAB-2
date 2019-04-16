%Joshua Enxing
%Tufts Univeristy
%MA226

%Gives approximate solution vector for u in a certain integral equation
%Input n is the number of subintervals used for the composite
%Simpson's rule
function u = solve_for_u(n)

x = linspace(0,1,n+1);
b = ((1/3)*((x.^2+1).^1.5-x.^3))';

A = zeros(n+1);
for i=1:n+1
    for j=1:n+1
            A(i,j) = (mod(j-1,2)*2+2)*(1/(3*n))*(x(i)^2+x(j)^2)^0.5;
    end
end

A(:,1) = A(:,1)*0.5;
A(:,n+1) = A(:,n+1)*0.5;

disp("Condition number of A for n = " + n + ": " + cond(A));
u = A\b;
plot(x,u,x,x);
legend('Computed Solution','True Solution');
title("N = " + n);
