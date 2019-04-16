%Joshua Enxing
%University of Connecticut
%MA5520
%Due 2/19/16

function [u_h, theta, h] = pw_linear_FEM (z, alpha, beta, a, b, f, m)

%pw_linear_FEM returns the values of u on the nodes given in the 
%vector z, with u(z_1) = alpha, u(z_(n+2)) = beta, calculated using
%piece-wise linear finite elements. 

%a, b, and f are the functions a(x), b(x), and f(x) in the differential 
%equation, and m is the degree of Gauss-Legendre quadrature desired for the 
%numerical integrations

n = length(z)-2;
h = zeros(n+1,1);

syms x;
a = sym(a);
b = sym(b);
f = sym(f);

for j=1:n+1
    h(j) = z(j+1) - z(j);
end

%generate the piece-wise linear finite elements
theta = generate_thetas(z,n);

%generate the matrix A as discussed in class
A = generate_matrix_A(a,b,z,h,n,m);
disp(A);

%generate the vector F as discussed in class, using the Gauss-Legendre
%weights and abcissas obtained from the gl_weight.m code provided in class
F = generate_vector_F(f,z,theta,alpha,beta,a,b,n,m);
disp(F);

%solve for the vector of values of u on our nodes in z
u_h = A\F;