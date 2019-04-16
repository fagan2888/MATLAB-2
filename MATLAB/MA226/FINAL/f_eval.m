function f = f_eval(p)

mat = load('vals.mat');
b = mat.b;
t = mat.t;

f = p(1)*exp(p(3).*t) + p(2)*exp(p(4).*t) - b;