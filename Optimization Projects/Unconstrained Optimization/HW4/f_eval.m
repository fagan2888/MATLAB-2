function fun = f_eval(A,b,x,m)

fun = zeros(m,1);

for i=1:m
    fun(i) = (1/10)*(0.5*(x')*A{i}*x - (b{i}')*x);
end





