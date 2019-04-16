t = (0 : .01 : 6)';

x1 = 0.5;
x2 = 0.5;
alpha1 = -0.3;
alpha2 = -0.4;

y = x1*exp(alpha1.*t) + x2*exp(alpha2.*t);

rng(8);
e = 10^(-4)*randn(601,1);

b = y + e;

save('vals','b','t');

p = [3 4 -1 -2]';
J = evaluate_jacobian(p);
c1 = cond(J);
p1 = lsqnonlin(@f_eval,p);
p1_err = [x1 x2 alpha1 alpha2]' - p1;

p =[3 4 -5 -6]';
J = evaluate_jacobian(p);
c2 = cond(J);
p2 = lsqnonlin(@f_eval,p);
p2_err = [x1 x2 alpha1 alpha2]' - p2;

p = [3 4 -0.3 -0.31]';
J = evaluate_jacobian(p);
c3 = cond(J);
p3 = lsqnonlin(@f_eval,p);
p3_err = [x1 x2 alpha1 alpha2]' - p3;
