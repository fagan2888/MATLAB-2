%solves w = L^{-1}r for input r that is a supervector of length 
%nK, where n = 100 and K = 40, and k = [k1 k2]
function y = get_y(k)

%Load and access data
[model] = prob_gen;
M = model.M;
A = model.A;
nx = model.nx;
nt = model.nt;
dt = model.dt;
y0 = model.y0;

y = zeros(nt*nx,1);

%Rename matrices
B = M/dt;
C = B + A(k);

%Find y_1 based on y0
y(1:nx) = C\(B*y0);

%Find y_i for the remaining time steps
for i=1:nt-1
    v = B*y((i-1)*nx+1:i*nx);
    y(i*nx+1:(i+1)*nx) = C\v;
end