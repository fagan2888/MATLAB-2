Jhat = @(k)calculate_jhat_k(k);

D = @(k,kt)calculate_direct_deriv(k,kt);

k = [0.1 0.1];

kt = [randn(1) randn(1)];

fdcheck(Jhat,D,k,kt);

fdcheck(Jhat,D,k,kt,4);