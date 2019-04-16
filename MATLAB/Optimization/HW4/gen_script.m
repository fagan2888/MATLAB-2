m = 10;
n = 2;
eta = 1;
[A,b] = generate_Ab(m,n,eta);

x = gradient_descent(A,b,m,[0 0]',2/101,100,1E-6);