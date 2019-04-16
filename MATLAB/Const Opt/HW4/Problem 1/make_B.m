%Makes B as in problem 2 for a given k
function B = make_B(k)

%Get y for this k
y = get_y(k);

%Load data
[model] = prob_gen;

%Access relevant data
nt = model.nt;
nx = model.nx;
A = model.A;

%Make and fill B as in problem 2
B = zeros(nx*nt,2);
for i=1:nt
    B((i-1)*nx+1:i*nx,1) = -A([1 0])*y((i-1)*nx+1:i*nx);
    B((i-1)*nx+1:i*nx,2) = -A([0 1])*y((i-1)*nx+1:i*nx);
end