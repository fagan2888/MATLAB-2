%Calculates the value of Jhat given a k value
function jhat = jhat(k)

%Load all necessary data
YHAT = load('yhat.mat');
YHAT = YHAT.yhat;

%Find the y_i values to use in Jhat for this particular k
Y = find_yk_yt_kt(k);

[model] = prob_gen;

M = model.M;
nt = model.nt;
dt = model.dt;

jhat = 0;

%Builds the sum for Jhat using the M inner product
for i=1:nt+1
    jhat = jhat + 0.5*dt*(Y(:,i)-YHAT(:,i))'*M*(Y(:,i)-YHAT(:,i));
end