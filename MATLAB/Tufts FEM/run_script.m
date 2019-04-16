close all;

%variable values
m =6;
h = 1/(2^m+2);
k = 2;
l = 2;
f = @(x,y) pi^2*(k^2+l^2)*sin(k*pi*x)*sin(l*pi*y);
tol = 1E-8;

x = linspace(0,1,1/h+1);
y = linspace(1,0,1/h+1);

F = make_vector(f,h);
%A = sparse(make_matrix(h));

t = ones(length(F),1);

P = cell(m,1);
A = cell(m,1);
for i=1:m
    P{m-i+1} = sparse(make_matrix_P((2^(m-i+1)+1)^2));
    h = 1/(2^(m-i+1)+2);
    if i==1
        A{m-i+1} = sparse(make_matrix(h));
    else
        A{m-i+1} = sparse((P{m-i+2})'*A{m-i+2}*(P{m-i+2}));
    end
end
h = 1/(2^m+2);

o = 1/h+1;
u1 = t;
v1 = zeros(o^2,1);
v2 = zeros(o^2,1);
counter = 0;
tic;
r = norm(F-A{m}*u1);
% while norm(F-A{m}*u1)>tol
%     u1 = vcycle(u1,F,A,P,1,1);
%     r_last = r;
%     r = norm(F-A{m}*u1);
%     counter = counter+1;
%     disp(r/r_last);
% end
u1 = FMG(F,A,P,1,1);
toc;

for i=2:o-1
        v1((i-1)*o+2:i*o-1) = u1((i-2)*(o-2)+1:(i-1)*(o-2));    
end

tic;
u2 = A{m}\F;
toc;

for i=2:o-1
    v2((i-1)*o+2:i*o-1) = u2((i-2)*(o-2)+1:(i-1)*(o-2));
end

w1 = v1;
w2 = v2;

%wrap vec into matrix
v1 = vec2mat(v1,1/h+1);
v2 = vec2mat(v2,1/h+1);

%plot solution
figure;
subplot(1,2,1);
surf(x,y,v1);
title('Multigrid');
subplot(1,2,2);
surf(x,y,v2);
title('Backslash');

e1 = 0;
e2 = 0;
[X,Y] = meshgrid(x,y);
X = X(:);
Y = Y(:);
DT = delaunayTriangulation(X,Y);
for i=1:length(w1)
    p = [DT.Points(i,1) DT.Points(i,2)];
    e1 = e1 + abs(w1(i)-sin(pi*k*p(1))*sin(pi*l*p(2)))^2;
    e2 = e2 + abs(w2(i)-sin(pi*k*p(1))*sin(pi*l*p(2)))^2;
end
e1 = num2str(sqrt(e1)*h);
e2 = num2str(sqrt(e2)*h);
disp("Multigrid error: " + e1 + ", Backslash error: " + e2 + ".");
disp("Vcycles: " + counter);
