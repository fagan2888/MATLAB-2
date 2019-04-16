%Outputs the value of the gradient of Jhat with respect to k
function D = gradient(k)

%Load and access data
[model] = prob_gen;
nt = model.nt;
nx = model.nx;

%Get y for this k
y = get_y(k);

%Make q as given in problem 3
q = make_q(k);

%Make B as in problem 2
B = make_B(k);

%Solve the adjoint equation with this q
lambda = action_inv_adjoint(q,k);

D = 0;

%Loop over all time steps
for i=1:nt
    D = D + lambda((i-1)*nx+1:i*nx)'*B((i-1)*nx+1:i*nx,:);
end
D = D';