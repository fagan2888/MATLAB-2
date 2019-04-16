%Performs Backward-Central finite difference method for u_t = -u_x. Input arguments are
%number of time nodes t_nodes, number of spatial nodes x_nodes, anonymous
%function initial condition init_cond (e.g. @(x)x.^2), anonymous function 
%left boundary condition left_bcs (e.g. @(t)t.^2), anonymous function right boundary
%condition right_bcs (e.g. @(t)t.^2), and number of time steps t_steps.

%Returns approximating solution vector u at time t_steps/(t_nodes-1).
function u = backward_central(t_nodes,x_nodes,init_cond,left_bcs,right_bcs,t_steps)

%Compute h_t, h_x, and lambda
h_t = 1/(t_nodes-1);
h_x = 1/(x_nodes-1);
lambda = h_t/(2*h_x);

%Evaluate the boundary condition functions at all time steps and put
%these values into vectors
[l_bcs,r_bcs] = make_bc_vecs(left_bcs,right_bcs,t_nodes,h_t);

%Evaluate the initial condition function for all x nodes and put this 
%into a vector
inits = make_init_cond_vec(init_cond,x_nodes,h_x);

%Make the P matrix
P = full(gallery('tridiag',x_nodes-2,-1*lambda,1,lambda));

%Initialize vector v of interior nodes
v = inits(2:length(inits)-1)';

%Initialize vector G 
G = zeros(length(v),1);

%Initialize solution vector u
u = zeros(length(v)+2,1);

%Make vector of x coordinate values
X = linspace(0,1,x_nodes);

%Iterate the crank-nicolson scheme, at each time step plotting u vs x
for i=1:t_steps
    G(1) = lambda*l_bcs(i+1);
    G(length(G)) = -1*lambda*r_bcs(i+1);
    
    v = v + G;
    v = P\v;
    
    u(1) = l_bcs(i+1);
    u(2:length(u)-1) = v;
    u(length(u)) = r_bcs(i+1);
    plot(X,u);
    pause(1)
end