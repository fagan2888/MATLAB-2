%Makes interpolation matrix P with n_rows rows and 1 + (n_rows-1)/2 columns
function P = make_matrix_P(n_rows)

n = n_rows;
l = log(n-1)/log(2);

P = zeros(n,2^(l-1)+1);

P(1,1) = 1;
P(2,1) = 1/2;
for i=2:2^(l-1)
    P(2*(i-1),i) = 1/2;
    P(2*(i-1)+1,i) = 1;
    P(2*(i-1)+2,i) = 1/2;
end
P(n-1,2^(l-1)+1) = 1/2;
P(n,2^(l-1)+1) = 1;