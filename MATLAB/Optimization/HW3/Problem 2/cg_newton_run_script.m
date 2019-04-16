for i=2:6
    n = 2^i;
    x = ones(n,1);
    x = x*(-1);
    z = newton_CG(x,0.1,1,0.5,0.9,100,1E-6,0);
end

n = 2^i;
x = ones(n,1);
x = x*(-1);
z = newton_CG(x,0.5,1,0.5,0.9,100,1E-6,1);
z = newton_CG(x,0.5,1,0.5,0.9,100,1E-6,2);