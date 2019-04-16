function f = f_eval(w0,w1,x0,x1)

f = zeros(4,1);
f(1) = w0 + w1 - 2;
f(2) = w0*x0 + w1*x1 - 2;
f(3) = w0*x0^2 + w1*x1^2 - (8/3);
f(4) = w0*x0^3 + w1*x1^3 - 4;