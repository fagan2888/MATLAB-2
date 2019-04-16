%Evaluate the initial condition function at all values of x
%and store these values in a vector
function inits = make_init_cond_vec(init_cond,x_nodes,h_x)

inits = linspace(0,1,x_nodes);

for i=1:x_nodes
    inits(i) = init_cond((i-1)*h_x);
end