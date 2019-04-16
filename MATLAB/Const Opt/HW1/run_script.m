[model] = prob_gen;

M = model.M;
dt = model.dt;
A = model.A;
t = model.t;
y0 = model.y0;
dx = model.dx;

y = y0;

k = [0.1 0.002];
B = M + dt*A(k);
for i = 1:20
    y = B\(M*y);
end

x = 0:0.01:1;
y = y';
plot(x,[y(100) y]);
    