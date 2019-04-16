function x_h = FMG(b,A_cell,P_cell,n_1,n_2)
n = length(b);
b_h = b;
r = log2(sqrt(n)-1);

if n == 25
    x_h = vcycle(zeros(25,1),b_h,A_cell,P_cell,n_1,n_2);
else
    P = P_cell{r};
    b_2h = (P')*b;
    %A_2h = (P')*A*P;
    x_2h = FMG(b_2h,A_cell,P_cell,n_1,n_2);
    x_h = P*x_2h;
end

x_h = vcycle(x_h,b_h,A_cell,P_cell,n_1,n_2);
x_h = vcycle(x_h,b_h,A_cell,P_cell,n_1,n_2);