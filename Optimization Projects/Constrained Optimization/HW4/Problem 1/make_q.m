%Makes the vector q for input k
function q = make_q(k)

%Load and access data
[model] = prob_gen;
M = model.M;
dt = model.dt;
nx = model.nx;
nt = model.nt;
yhat = load('yhat.mat');
yhat = yhat.yhat;

%Find y for this k
y = get_y(k);

%Make and fill q
q = zeros(nx*nt,1);
for i=1:nt
    q((i-1)*nx+1:i*nx) = dt*M*(y((i-1)*nx+1:i*nx)-yhat(:,i+1));
end

