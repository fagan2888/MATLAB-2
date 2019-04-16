%solves s = L^{-*}q for input q that is a supervector of length 
%nK, where n = 100 and K = 40, and k = [k1 k2]
function s = action_inv_adjoint(q,k)

%Load and access data
[model] = prob_gen;
M = model.M;
A = model.A;
nx = model.nx;
nt = model.nt;
dt = model.dt;

s = zeros(nt*nx,1);

%Rename matrices
B = M/dt;
C = B + A(k);

%Find s for the Kth time step
s((nt-1)*100+1:nt*100) = C\q((nt-1)*100+1:nt*100);

%Find s for remaining time steps (in backward time)
for i=2:nt
    j = nt-i+1;
    v = B*s(j*100+1:j*100+100) + q((j-1)*100+1:j*100);
    s((j-1)*100+1:(j-1)*100+100) = C\v;
end


