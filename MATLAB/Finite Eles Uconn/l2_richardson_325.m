%Joshua Enxing
%University of Connecticut
%MA5520
%Due 3/25/16

function [l2_err_12,l2_err_24,z2,z4,r1,r2,r4] = l2_richardson_325(h,alpha,beta,a,b,f,m)

d = h/4;
syms x;
z4 = linspace(0,1,round((1/d)+1));
[u_4h,theta_4h,h] = pw_linear_FEM_325(z4,alpha,beta,a,b,f,m);

z2 = linspace(0,1,round((1/(2*d))+1));
[u_2h,theta_2h,h2] = pw_linear_FEM_325(z2,alpha,beta,a,b,f,m);

z1 = linspace(0,1,round((1/(4*d))+1));
[u_h,theta_h,h] = pw_linear_FEM_325(z1,alpha,beta,a,b,f,m);

r1 = zeros(length(u_h)+2,1);
r2 = zeros(length(u_2h)+2,1);
r4 = zeros(length(u_4h)+2,1);

r1(1) = alpha;
r2(1) = alpha;
r4(1) = alpha;

r1(length(u_h)+2) = beta;
r2(length(u_2h)+2) = beta;
r4(length(u_4h)+2) = beta;

r1(2:length(u_h)+1) = u_h;
r2(2:length(u_2h)+1) = u_2h;
r4(2:length(u_4h)+1) = u_4h;

l2_err_12 = 0;
l2_err_24 = 0;

for i=1:length(z2)-1
    c = ceil(i/2);
    [y,w] = gl_weight(z2(i),z2(i+1),m);
    for j=1:length(y)
        l2_err_24 = l2_err_24 + subs(((r1(c)*theta_h(c)+r1(c+1)*theta_h(c+1))-(r2(i)*theta_2h(i)+r2(i+1)*theta_2h(i+1)))^2,x,y(j))*w(j);
    end
end

for i=1:length(z4)-1
    c = ceil(i/2);
    [y,w] = gl_weight(z4(i),z4(i+1),m);
    for j=1:length(y)
        l2_err_12 = l2_err_12 + subs(((r2(c)*theta_2h(c)+r2(c+1)*theta_2h(c+1))-(r4(i)*theta_4h(i)+r4(i+1)*theta_4h(i+1)))^2,x,y(j))*w(j);
    end
end

l2_err_12 = double(sqrt(l2_err_12));
l2_err_24 = double(sqrt(l2_err_24));
