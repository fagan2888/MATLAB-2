function [y,u,D] = chebyshev_spectral_dirichlet(init,N,t_f,dt)

syms x;
y = make_collocation_points(N);
D = make_diff_matrix(y);
D2 = D^2;
u = zeros(N,1);
num_steps = t_f/dt;

for i=1:N
    u(i) = init(y(i));
end

for i=1:num_steps
    k1 = D*u;
    u_k2 = u + 0.5*dt*k1;
    k2 = D*u_k2;
    u_k3 = u + 0.5*dt*k2;
    k3 = D*u_k3;
    u_k4 = u + dt*k3;
    k4 = D*u_k4;
    
    du = k1+2*k2+2*k3+k4;
    
    u = u + (dt/6)*du;
    u(N) = bdyl(dt*i);
end

plot(y,u);
