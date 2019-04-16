function norm = L2_norm_function(f,a,b)

f = @(x)f(x).*f(x);
norm = sqrt(integral(f,a,b));