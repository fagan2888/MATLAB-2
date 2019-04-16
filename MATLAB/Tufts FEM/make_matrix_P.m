%Makes interpolation matrix P with n_rows rows and 1 + (n_rows-1)/2 columns
function P = make_matrix_P(n_rows)

n = n_rows;
if mod(sqrt(n),2)~=0
    m = ((sqrt(n)+1)/2)^2;
else
    m = (sqrt(n)/2)^2;
end

if n==9 
    P = [1/2 1/2 0 1/2 1 1/2 0 1/2 1/2]';
else
    P = zeros(n,m);
    
    a = sqrt(n);
    b = sqrt(m);
    
    A = zeros(a,b);
    B = zeros(a,b);
    D = zeros(a,b);
    
    for i=1:b-1
        A(1+(i-1)*2,i) = 1;
        A(2+(i-1)*2,i:i+1) = [1/2 1/2];
        
        B(1+(i-1)*2:2+(i-1)*2,i) = [1/2 1/2]';
    end
    A(a,b) = 1;
    B(a,b) = 1/2;
    
    C = rot90(B);
    C = rot90(C);
    
    for i=1:b-1
        P(1+2*a*(i-1):2*a*i,1+b* (i-1):b*(i+1)) = [A D;B C];
    end
    
    P(n-a+1:n,m-b+1:m) = A;
    
end

