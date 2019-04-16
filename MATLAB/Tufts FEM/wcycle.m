%Does one iteration of wcycle
function x = wcycle(x_h,b_h,A_cell,P_cell,nu_1,nu_2)

x = x_h;
b = b_h;
n = length(x);
r = log2(sqrt(n)-1);
A = A_cell{r};

%Relax nu_1 times with gauss-seidel
for i=1:nu_1
    x = gauss_seidel(A,x,b);
    %x = (b-A*x+diag(diag(A))*x)./diag(A);
end

%make interpolation matrix P
P = P_cell{r};

%Restrict A
%A_2h = (P')*A*P;
A_2h = A_cell{r-1};

%Restrict residual
res = (P')*(b - A*x);

%If on next-to-coarsest grid, restrict and direct-solve, else
%recursively call vcycle solving Ae = r on the coarser grid
if n==25
    e = A_2h\res;
else 
    temp = zeros(length(res),1);
    e = vcycle(temp,res,A_cell,P_cell,nu_1,nu_2);
    e = vcycle(e,res,A_cell,P_cell,nu_1,nu_2);
end

%Correct
x = x + P*e;

%Relax nu_2 times with gauss-seidel from hw5
for i=1:nu_2
    x = gauss_seidel(A,x,b);
    %x = (b-A*x+diag(diag(A))*x)./diag(A);
end