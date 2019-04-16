%Joshua Enxing
%University of Connecticut
%MA5520
%Due 2/19/16

%generate the l2 error between the true solution u and the FEM solution
function [err,r] = l2_error(u,alpha, beta, u_h,theta,v,m)

syms x;

%make r, vector of values of FEM solution on all nodes. so we are just
%adding in the BCs here
r = zeros(length(u_h)+2,1);
r(1) = alpha;
r(length(u_h)+2) = beta;
r(2:length(u_h)+1) = u_h;

err = 0;

%integrate (u-u_h)^2 on each element and sum all of these values
for i=1:length(v)-1
    [y,w] = gl_weight(v(i),v(i+1),m);
    for j=1:length(y)
        err = err + subs((u-(r(i)*theta(i)+r(i+1)*theta(i+1)))^2,x,y(j))*w(j);
    end
end

err = double(sqrt(err));
