K = zeros(20,1);
iters = zeros(20,1);

for i = 1:20
    v = rand(2,1);
    v(1) = 2*i*v(1);
    A = diag(v);
    K(i) = cond(A);
    [x,iters(i)] = steep_desc(A,[1 1]',1,[0 0]',1E-6,10000000);
end

scatter(K,iters);
