%Does one iteration of vcycle
function x = vcycle(x_h,b_h,A_cell,P_cell,nu_1,nu_2)

x = x_h;
b = b_h;
n = length(x);
r = log2(sqrt(n)-1);
A = A_cell{r};

%Relax nu_1 times with Jacobi
for i=1:nu_1
    %x = gauss_seidel(A,x,b);
    x = x + tril(A)\(b-A*x);
    %x = (2/3)*(b-A*x+diag(diag(A))*x)./diag(A)+(1/3)*x;
end

%make interpolation matrix P
P = P_cell{r};

%If on next-to-coarsest grid, restrict and direct-solve, else
%recursively call vcycle solving Ae = r on the coarser grid
if n==25
    P = P_cell{2};
    
    %Restrict A
    %A_2h = (P')*A*P;
    A_2h = A_cell{1};

    %Restrict residual
    res = (P')*(b - A*x);
    
    e = A_2h\res;
else 
    %Restrict A
    %A_2h = (P')*A*P;

    %Restrict residual
    res = (P')*(b - A*x);
    
    temp = zeros(length(res),1);
    e = vcycle(temp,res,A_cell,P_cell,nu_1,nu_2);
end

%Correct
x = x + P*e;

%Relax nu_2 times with GS
for i=1:nu_2
    %x = gauss_seidel(A,x,b);
    x = x + tril(A)\(b-A*x);
    %x = (2/3)*(b-A*x+diag(diag(A))*x)./diag(A)+(1/3)*x;
end