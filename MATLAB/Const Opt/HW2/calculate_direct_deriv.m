%Calculates the value of the directional derivative of Jhat with respect to
%k in the direction of kt
function deriv = calculate_direct_deriv(k,kt)

%Load necessary data
YHAT = load('yhat.mat');
YHAT = YHAT.yhat;

%Find all y_i's and {y tilde}_i's for all time steps
[Y,YT] = find_yk_yt_kt(k,kt);

[model] = prob_gen;

M = model.M;
nt = model.nt;
dt = model.dt;

deriv = 0;

%Build directional derivative using M inner product
for i=1:nt+1
    deriv = deriv + dt*(Y(:,i)-YHAT(:,i))'*M*(YT(:,i));
end