%solves w = L^{-1}r for input r that is a supervector of length 
%nK, where n = 100 and K = 40, and k = [k1 k2]
function w = action_L_inv(r,k)

%Load and access data
[model] = prob_gen;
M = model.M;
A = model.A;
nx = model.nx;
nt = model.nt;
dt = model.dt;

w = zeros(nt*nx,1);

%Renam matrices
B = M/dt;
C = B + A(k);

%Find w for the first time step
w(1:nx) = C\r(1:nx);

%Find w for the remaining time steps
for i=1:nt-1
    v = B*w((i-1)*nx+1:i*nx) + r(i*nx+1:(i+1)*nx);
    w(i*nx+1:(i+1)*nx) = C\v;
end