m = 10;
n = 2;
eta = 1;
[A,b] = generate_Ab(m,n,eta);

M = zeros(2);
v = zeros(2,1);

for i=1:10
    M = M + A{i};
    v = v + b{i};
end

sol = M\v;
disp(sol);

[x,err1] = gradient_descent(A,b,m,[0 0]',2/101,500,sol,1E-6);

figure;
plot(err1);
title('GD');

err1 = log(err1);
err1_y = err1(2:length(err1));

figure;
plot(err1(1:length(err1)-1),err1_y);
title('GD error');

[x,err2] = inc_gradient_descent(A,b,m,[0 0]',2/101,50000,sol,1E-6);

l = err2;
figure;
plot(err2);
title('IGD');

err2 = log(err2);
err2_y = err2(2:length(err2));

figure;
plot(err2(1:length(err2)-1),err2_y);
title('IGD error');

[x,err3] = stoch_gradient_descent(A,b,m,[0 0]',2/101,5000,sol,1E-6);

figure;
plot(err3);
title('SGD');

err3 = log(err3);
err3_y = err3(2:length(err3));

figure;
plot(err3(1:length(err3)-1),err3_y);
title('SGD error');