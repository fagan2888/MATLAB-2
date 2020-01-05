%
% Centered Finite Difference Derivative Checker
%
% Computes the "Finite Difference V-curve" by the 2nd order and 4th order
% finite difference approximation to the derivative. Plots the
% error between the provided directional derivative, fp, and the finite
% difference approximation of provided function, f. If fp is the correct
% directional derivative of f the error will form the "Finite Difference V"
%
% For a lookup table of finite difference coefficients, see:
% https://en.wikipedia.org/wiki/Finite_difference_coefficient
%
%
% input
% f anonymous function (scalar-, vector-, or matrix-valued)
% fp anonymous function, directional derivative of f, with two
% arguments, fp(x, xp) = D_x[f(x)](xp)
% (scalar-, vector-, or matrix-valued)
% x center point for the derivative check
% xp direction of perturbation of x
% ord (optional, default = 2) order of the approximation, ord = {2,4}
%
%
% output: min_error minimum FD error in the "Finite Difference V"

function [ min_error ] = fdcheck( f, fp, x, xp, ord )

%Default ord is 2
if nargin < 5
    ord = 2;
end

%Make grid of epsilons
grid = logspace(-8,-2,100);

%Make error vector
err_vect = zeros(100,1);

%When ord = 2, compute the finite diff approximations and store the error
%with the actual in err_vect
if ord == 2
    for i=1:100
        approx = (f(x+grid(i)*xp)-f(x-grid(i)*xp))/(2*grid(i));
        err_vect(i) = abs(approx-fp(x,xp));
    end
end


%When ord = 4, compute the finite diff approximations and store the error
%with the actual in err_vect
if ord == 4
    for i=1:100
        approx = (-f(x+2*grid(i)*xp) + 8*f(x+grid(i)*xp) - 8*f(x-grid(i)*xp) + f(x-2*grid(i)*xp))/(12*grid(i));
        err_vect(i) = abs(approx-fp(x,xp));
    end
end

ord = num2str(ord);

%Plot abs error vs value of epsilon
figure;
loglog(grid,err_vect);
axis([10E-8 10E-2 0 inf]);
xlabel('Epsilon');
ylabel('Absolute error');
title("Error in FD Approximation, ord = "+ord);
min_error = min(err_vect);


