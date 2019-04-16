%Make a random (symmetric) Toeplitz plus Hankel matrix.
%n is dimension of matrix; vectors c, q, and p have random 
%integer values from 1 to m.

function [M,T,H] = RandomTHMatrix(n,m)

%Make random vector c
c = randi(m,n,1);
c = c';

%Make random vector r of length n-1, since the first entry of q must 
%agree with the last entry of c for H to be hankel
r = randi(m,n-1,1);
r = r';
q = zeros(1,n);
q(1) = c(n);
q(1,2:n) = r;

%Make random vector for the toeplitz matrix
p = randi(m,n,1);
p = p';

%Form hankel matrix
H = hankel(c,q);

%Form symmetric toeplitz matrix
T = toeplitz(p);

M = T+H;