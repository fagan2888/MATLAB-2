%Joshua Enxing
%University of Connecticut
%MA5520
%Due 3/25/16

%generates the vector F as discussed in class, but with only the first and
%last entries being nonzero, due to zero forcing in our new problem
function F = generate_vector_F_325(f, z, theta, alpha, beta, a, b, n, m)
syms x;
F = zeros(n,1);

%generate weights and abcissas for extra term of F_1 and F_N
[y1,w1] = gl_weight(z(1),z(2),m);
[yn,wn] = gl_weight(z(n+1),z(n+2),m);

for j=1:length(y1)
    %integrate -alpha*[a*(theta_0)'*(theta_1)' + b*(theta_0)'*theta_1] on
    %[x_0,x_1] using g-l quadrature weights and abcissas
    F(1) = F(1) - alpha*subs(a*(-1/(z(2))^2)+(b*(x-z(1)))/(-(z(1)-z(2))^2),x,y1(j))*w1(j);
    
    %integrate -alpha*[a*(theta_n+2)'*(theta_n+1)' + b*(theta_n+2)'*theta_n+1] on
    %[x_0,x_1] using g-l quadrature weights and abcissas
    F(n) = F(n) - beta*subs(a*(-1/(1-z(n+1))^2)+(b*(x-z(n+2)))/(-(z(n+1)-z(n+2))^2),x,yn(j))*wn(j);
end