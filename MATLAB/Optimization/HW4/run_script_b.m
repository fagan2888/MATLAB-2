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

[x,err] = reshuffling_gradient_descent(A,b,m,[0 0]',2/101,5000,sol,1E-6);

figure;
plot(err);
title('SGD Random Reshuffling');

err = log(err);
err_y = err(2:length(err));

figure;
plot(err(1:length(err)-1),err_y);
title('SGD Random Reshuffling Error');