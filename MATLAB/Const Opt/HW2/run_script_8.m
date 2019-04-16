%Make function for Jhat
Jhat = @(k)calculate_jhat_k(k);

%Make function for directional derivative
D = @(k,kt)calculate_direct_deriv(k,kt);

%Fix k
k = [0.1 0.1];

%Find random k tilde direction
kt = [randn(1) randn(1)];

%Check second-order FD
fdcheck(Jhat,D,k,kt);

%Check fourth-order FD
fdcheck(Jhat,D,k,kt,4);