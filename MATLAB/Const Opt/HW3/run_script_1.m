%Make random r and q
r = rand(4000,1);
q = rand(4000,1);

%Make arbitrary k
k = rand(2,1);

%Find L^{-*}q for this k
s = action_inv_adjoint(q,k);

%Find L^{-1}r for this k
w = action_L_inv(r,k);

format long

%Display error
disp("Abs value of <L^{-1}r,q> - <r,L^{-*}q>: " + abs(w'*q-r'*s));