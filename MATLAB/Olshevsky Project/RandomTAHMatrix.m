%Make a random (symmetric) Toeplitz plus Hankel matrix.
%n is dimension of matrix; vectors c, q, and p have random 
%integer values from 1 to m.

function [M,T,H] = RandomTAHMatrix(n,m)

%Make random vector c
c = randi(m,n,1);
c = c';

%Make vector q such that q(i) = c(n-(i-1)) for i from 1 to n so Hankel 
%is antisymmetric
q = zeros(1,n);
for i=1:n
    q(i) = c(n - (i-1));
end

%Make random vector for the toeplitz matrix
p = randi(m,n,1);
p = p';

%Form antisymmetric hankel matrix
H = hankel(c,q);

%Form symmetric toeplitz matrix
T = toeplitz(p);

M = T+H;