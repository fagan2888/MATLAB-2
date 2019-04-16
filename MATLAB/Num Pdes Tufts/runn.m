function [error] = runn(u_true,rhs)

syms x;
sp = sym(zeros(4,1));
fe = sym(zeros(4,1));
error = zeros(5,4);
for i=1:5
    sp(i) = spectral_legendre(2*2^(i-1),rhs);
    fe(i) = FE_neumann(2*2^(i-1),rhs);
    error(i,1) = double(int((u_true-sp(i))^2,x,-1,1));
    error(i,2) = double(int((u_true-fe(i))^2,x,-1,1));
    error(i,3) = double(int((diff(u_true-sp(i)))^2+(u_true-sp(i))^2,x,-1,1));
    error(i,4) = double(int((diff(u_true-fe(i)))^2+(u_true-fe(i))^2,x,-1,1));
end
