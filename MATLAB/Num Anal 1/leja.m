%Joshua Enxing%
%University of Connecticut%
%MA5510%

function [x,f]=leja(x,f)

n = length(x);

vx = zeros(1,n);
vf = zeros(1,n);
m = zeros(1,n);

%Find first value in leja ordering% 
[M,I] = max(abs(x));
m(1) = I;
vx(1) = x(m(1));
vf(1) = f(m(1));

%Find ordering for the rest of the x's%
for k=1:n-1
    M = (-1)*eye(k+1); %M is subtraction matrix%
    M(k+1,k+1)=1;
    R = zeros(k+1,n);%R is matrix of columns of previous x_k's with x_i at bottom%
                     %only bottom row of R differs from column to column%
    for j=1:k
        M(j,k+1)=1;
        R(j,:)=x(m(j));
    end
    for i=1:n
        R(k+1,i)=x(i);
    end
    
    %Matrix of differences, columns specific to each x(i)%
    D = M*R;
    D(k+1,:)=1;
    D = abs(D);
    
    %B is product of differences between previous x_k's and current x_i% 
    B = prod(D);
    [q,I]= max(B);
    m(k+1)= I;
    vx(k+1) = x(m(k+1));
    vf(k+1) = f(m(k+1));  
end

x = vx';
f = vf';



