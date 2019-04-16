%Create function handle for objective function
f = @(k) calculate_jhat_k(k);

%Create function handle for directional derivative
fp = @(k,kt) dir_deriv(k,kt);

%Make random vectors k and kt
k = rand(2,1);
kt = rand(2,1);

%Run finite difference checker
e = fdcheck(f,fp,k,kt);