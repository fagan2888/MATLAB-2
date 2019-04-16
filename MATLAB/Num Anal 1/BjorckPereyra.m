%Joshua Enxing%
%University of Connecticut%
%MA5510%

function a = BjorckPereyra(x,f)

a = f;

%Create the n-1 right L matrices%
n = length(x);
LR = cell(n-1);

for i=1:n-1
    LR{i} = eye(n);
    for k=i+1:n
        LR{i}(k,i) = -1;
    end
end

%create difference vectors%
v = cell(n-1);

for k = 1:n-1
    v{k}=LR{k}*x;
end


%Create n-1 left L matrices%
LL=cell(n-1);
for i=1:n-1
    M = eye(n);
    for k=i+1:n
        M(k,k) = 1/(v{i}(k));
    end
    LL{i} = M;
end

%Create n-1 U matrices%
N = (-1)*x;
U=cell(n-1);
for i=1:n-1
    U{i} = eye(n);
    for k=i:n-1
        U{i}(k,k+1)= N(i);
    end
end

for i=1:n-1
    a = LL{i}*LR{i}*a;
end

for i=1:n-1
    a = U{n-i}*a;
end

