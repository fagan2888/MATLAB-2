function alpha = compute_alpha(A,b,s,x)

alpha = ((s')*b-(s')*(A')*A*x)/((s')*(A')*A*s);