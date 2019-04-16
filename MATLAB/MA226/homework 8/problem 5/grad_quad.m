function grad = grad_quad(A,b,x)

grad = (A')*(A*x) - b;
