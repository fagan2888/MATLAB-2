function err = err_for_eps(e,a,b)

[model] = prob_gen;

M = model.M;
dt = model.dt;
A = model.A;
t = model.t;
y0 = model.y0;
dx = model.dx;

y = y0;
Y = zeros(100,20);
k = [0.1 0.002];
kt = [a b];

B = M + dt*A(k);
for i = 1:20
    y = B\(M*y);
    Y(:,i) = y;
end

yt = zeros(100,1);
YT = zeros(100,20);
for i = 1:20
    C = (M*yt - dt*A(kt)*Y(:,i));
    yt = B\C;
    YT(:,i) = yt;
end

y = y0;
Y_plus = zeros(100,20);
B = M + dt*A(k+e*kt);
for i = 1:20
    y = B\(M*y);
    Y_plus(:,i) = y;
end

y = y0;
Y_minus = zeros(100,20);
B = M + dt*A(k-e*kt);
for i = 1:20
    y = B\(M*y);
    Y_minus(:,i) = y;
end

approx = (1/(2*e))*(Y_plus-Y_minus);

err = norm(YT-approx);

disp("Norm of diff for e = " + e + " : " + err);


