%Joshua Enxing
%University of Connecticut
%MA5520
%Due 2/19/16

%generates the matrix A as discussed in class
function A = generate_matrix_A(a,b,z,h,n,m)

syms x;

%K, stiffness matrix
c = zeros(1,n);
d = zeros(1,n-1);

%M, mass matrix
f = zeros(1,n);
g = zeros(1,n-1);

%Generate the values of the diagonals and off-diagonals of the matrices
%M and K using quadrature for the integrals
for i=1:n+1
    [y,w] = gl_weight(z(i),z(i+1),m);
    if i==1
        for j=1:length(y)
            c(i) = c(i) + subs(a/((z(i+1)-z(i))^2),x,y(j))*w(j);
            f(i) = f(i) + subs(b*((x-z(i))/(z(i+1)-z(i)))^2,x,y(j))*w(j);
        end
        
    elseif i==n+1
        for j=1:length(y)
            c(i-1) = c(i-1) + subs(a/((z(i+1)-z(i))^2),x,y(j))*w(j);
            f(i-1) = f(i-1) + subs(b*((x-z(i))/(z(i+1)-z(i)))^2,x,y(j))*w(j);
        end
    
    else
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
    
            
%This following commented code was using the approximations mentioned in the slides, 
%which did not produce error on the order of roundoff for part (b).

%for i = 1:n
%    if i == 1
%        c(1) = (subs(a,x,(z(2)-0.5*h(1)))/h(1)) + (subs(a,x,(z(2)+0.5*h(2)))/h(2));
%        f(1) = (subs(b,x,(z(2)-0.5*h(1)))*h(1))/3 + (subs(b,x,(z(2)+0.5*h(2)))*h(2))/3;
%    
%    else
%        c(i) = (subs(a,x,(z(i+1)-0.5*h(i)))/h(i)) + (subs(a,x,(z(i+1)+0.5*h(i+1)))/h(i+1));
%        f(i) = (subs(b,x,(z(i+1)-0.5*h(i)))*h(i))/3 + (subs(b,x,(z(i+1)+0.5*h(i+1)))*h(i+1))/3;
%        d(i-1) = -subs(a,x,(z(i+1)-0.5*h(i)))/h(i);
%        g(i-1) = (subs(b,x,(z(i+1)-0.5*h(i)))*h(i))/6;
%    end   
%end

K = gallery('tridiag',d,c,d);
M = gallery('tridiag',g,f,g);

A = M+K;
