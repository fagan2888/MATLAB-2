function fun = f_eval(x)

[ A, b ] = generate_Ab(10, 2, 1);
fun = zeros(10,1);

for i=1:10
    fun(i) = (1/10)*(0.5*(x')*A{i}*x - (b{i}')*x);
end