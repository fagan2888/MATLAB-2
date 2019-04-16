%Calculates matrices Y and YT whose columns are the values of y and y tilde
%for each time step, based on the particular inputs for k and k tilde
function [Y,YT] = find_yk_yt_kt(k,kt)

%Load all necessary data
[model] = prob_gen;

M = model.M;
dt = model.dt;
A = model.A;
nt = model.nt;
y0 = model.y0;

y = y0;

Y = zeros(100,nt+1);

%Recycled code from hw1 since system is the same
B = M + dt*A(k);
Y(:,1) = y;
for i = 2:nt+1
    y = B\(M*y);
    Y(:,i) = y;
end

%If leave k tilde argument empty (i.e. only using this to calculate Jhat,
%don't do this step to waste time/computation
if nargin == 2
    %Recycled code from hw1 since system is the same
    YT = zeros(100,nt+1);
    B = M + dt*A(k);
    for i = 2:nt+1
        C = (M*YT(:,i-1) - dt*A(kt)*Y(:,i));
        YT(:,i) = B\C;
    end
end