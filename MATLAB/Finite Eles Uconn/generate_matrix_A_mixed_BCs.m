%Joshua Enxing
%University of Connecticut
%MA5520
%Due 2/19/16

%generates the matrix A as discussed in class
function A = generate_matrix_A_mixed_BCs(a,b,z,h,n,m)

syms x;

%K, stiffness matrix
c = zeros(1,n);
d = zeros(1,n-1);

%M, mass matrix
f = zeros(1,n);
g = zeros(1,n-1);

%Generate the values of the diagonals and off-diagonals of the matrices
%M and K using quadrature for the integrals
for i=1:n
    if i==1
        [y,w] = gl_weight(z(i),z(i+1),m);
        for j=1:length(y)
            c(i) = c(i) + subs(a/((z(i+1)-z(i))^2),x,y(j))*w(j);
            f(i) = f(i) + subs(b*((x-z(i))/(z(i+1)-z(i)))^2,x,y(j))*w(j);
        end
    
    else
        [y,w] = gl_weight(z(i),z(i+1),m);
        for j=1:length(y)
            c(i) = c(i) + subs(a/((z(i+1)-z(i))^2),x,y(j))*w(j);
            c(i-1) = c(i-1) + subs(a/((z(i+1)-z(i))^2),x,y(j))*w(j);
            d(i-1) = d(i-1) + subs(a/(-(z(i+1)-z(i))^2),x,y(j))*w(j);
            f(i) = f(i) + subs(b*((x-z(i))/(z(i+1)-z(i)))^2,x,y(j))*w(j);
            f(i-1) = f(i-1) + subs(b*((x-z(i))/(z(i+1)-z(i)))^2,x,y(j))*w(j);
            g(i-1) = g(i-1) + subs((b*(x-z(i+1))*(x-z(i)))/(-(z(i+1)-z(i))^2),x,y(j))*w(j);
        end
    end
    
end

K = gallery('tridiag',d,c,d);
M = gallery('tridiag',g,f,g);

A = M+K;