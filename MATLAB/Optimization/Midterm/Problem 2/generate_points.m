function [x,y] = generate_points

rand('seed',314);
x = linspace(0,1,30);
y = 2*x.^2-3.*x+1+0.05*rand(size(x)) ; 