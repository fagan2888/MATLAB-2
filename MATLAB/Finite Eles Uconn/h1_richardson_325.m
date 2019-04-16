%Joshua Enxing
%University of Connecticut
%MA5520
%Due 3/25/16

function [l2_p_325,h1_p_325,h] = h1_richardson_325(h,alpha,beta,a,b,f,m)

syms x;

[l2_err_12,l2_err_24,z2,z4,r1,r2,r4] = l2_richardson_325(h,alpha,beta,a,b,f,m);

l2_p_325 = log(l2_err_24/l2_err_12)/log(2);

h1_err_12 = 0;
h1_err_24 = 0;

for i=1:length(z2)-1
    c = ceil(i/2);
    d = z2(i+1)-z2(i);
    h1_err_24 = h1_err_24 + d*(((r1(c+1)-r1(c))/(2*d))-((r2(i+1)-r2(i))/d))^2;
end

for i=1:length(z4)-1
    c = ceil(i/2);
    d = z4(i+1)-z4(i);
    h1_err_12 = h1_err_12 + d*(((r2(c+1)-r2(c))/(2*d))-((r4(i+1)-r4(i))/d))^2;
end

h1_err_12 = h1_err_12 + (l2_err_12)^2;
h1_err_12 = double(sqrt(h1_err_12));

h1_err_24 = h1_err_24 + (l2_err_24)^2;
h1_err_24 = double(sqrt(h1_err_24));

h1_p_325 = log(h1_err_24/h1_err_12)/log(2);