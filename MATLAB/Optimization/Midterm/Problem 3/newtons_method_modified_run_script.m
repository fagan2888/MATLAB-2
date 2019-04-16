syms x y;
f = x.^4 + y.^2;

disp("Initial guess [0.5 0.5]'");
[z,e,iters] = newtons_method_modified(f,[x y],[0.5 0.5]',1000,1E-6);
e = log(e);
x = e(1:length(e)-1);
y = e(2:length(e));

figure
plot(x,y);
title("Initial Guess [0.5 0.5]'");
legend('log(e_i) vs log(e_{i+1})');
disp(" ");

syms x y;
f = x.^4 + y.^2;

disp("Initial guess [0.1 0.1]'");
[z,e,iters] = newtons_method_modified(f,[x y],[0.1 0.1]',1000,1E-6);
e = log(e);
x = e(1:length(e)-1);
y = e(2:length(e));

figure
plot(x,y);
title("Initial Guess [0.1 0.1]'");
legend('log(e_i) vs log(e_{i+1})');