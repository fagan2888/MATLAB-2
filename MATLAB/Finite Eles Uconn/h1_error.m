%Joshua Enxing
%University of Connecticut
%MA5520
%Due 2/19/16

%generate the h1 error between the true solution u and the FEM solution
function [h1err,l2err] = h1_error(u,alpha, beta, u_h,theta,h,v,m)

syms x;

%compute the l2 error between u and the FEM solution. here r is
%u_h but with two extra entries for the BCs
[l2err,r] = l2_error(u,alpha, beta, u_h,theta,v,m);
l = (l2err)^2;

%find the derivative of the true solution
u_prime = diff(u);

h1err = 0;

%integrate (u'-(u_h)')^2 on each element and sum all of these values
for i=1:length(v)-1
    [y,w] = gl_weight(v(i),v(i+1),m);
    for j=1:length(y)
        h1err = h1err + subs((u_prime-(r(i+1)/h(i)-r(i)/h(i)))^2,x,y(j))*w(j);
    end
end

%add square of l2 error to the integral of (u'-(u_h)')^2
h1err = l + h1err;

h1err = double(sqrt(h1err));


