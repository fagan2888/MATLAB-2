function f = compute_f(A,b,x,c)

f = 0.5*(x')*(A')*A*x + (x')*b + c;