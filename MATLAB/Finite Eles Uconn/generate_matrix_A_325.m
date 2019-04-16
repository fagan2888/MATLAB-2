%Joshua Enxing
%University of Connecticut
%MA5520
%Due 3/25/16

%generates the matrix A as discussed in class, except with the mass matrix
%being different for this problem (here it is M_ij =
%integral(b*(theta_i)'*(theta_j)) over the domain
function A = generate_matrix_A_325(a,b,z,h,n,m)

syms x;

%K, stiffness matrix
c = zeros(1,n);
d = zeros(1,n-1);

%M, mass matrix
f = zeros(1,n);
g = zeros(1,n-1);
s = zeros(1,n-1);

%Generate the values of the diagonals and off-diagonals of the matrices
%M and K using quadrature for the integrals
for i=1:n+1
    [y,w] = gl_weight(z(i),z(i+1),m);
    if i==1
        for j=1:length(y)
            c(i) = c(i) + subs(a/((z(i+1)-z(i))^2),x,y(j))*w(j);
            f(i) = f(i) + subs((-b*(x-z(i)))/(z(i+1)-z(i))^2,x,y(j))*w(j);
        end
        
    elseif i==n+1
        for j=1:length(y)
            c(i-1) = c(i-1) + subs(a/((z(i+1)-z(i))^2),x,y(j))*w(j);
            f(i-1) = f(i-1) + subs((b*(x-z(i)))/(z(i+1)-z(i))^2,x,y(j))*w(j);
        end
    
    else
        for j=1:length(y)
            c(i) = c(i) + subs(a/((z(i+1)-z(i))^2),x,y(j))*w(j);
            c(i-1) = c(i-1) + subs(a/((z(i+1)-z(i))^2),x,y(j))*w(j);
            d(i-1) = d(i-1) + subs(a/(-(z(i+1)-z(i))^2),x,y(j))*w(j);
            f(i) = f(i) + subs((-b*(x-z(i)))/(z(i+1)-z(i))^2,x,y(j))*w(j);
            f(i-1) = f(i-1) + subs((b*(x-z(i)))/(z(i+1)-z(i))^2,x,y(j))*w(j);
            g(i-1) = g(i-1) + subs((b*(x-z(i+1)))/(-(z(i+1)-z(i))^2),x,y(j))*w(j);
            s(i-1) = s(i-1) + subs((b*(x-z(i)))/(-(z(i+1)-z(i))^2),x,y(j))*w(j);
        end
    end
    
end

K = gallery('tridiag',d,c,d);
M = gallery('tridiag',s,f,g);

A = M+K;