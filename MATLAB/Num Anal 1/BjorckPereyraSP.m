%Joshua Enxing%
%University of Connecticut$
%MA5510%

function a = BjorckPereyraSP(x,f)

a = single(f);

%Create the n-1 right L matrices%
n = length(x);
LR = cell(n-1);

for i=1:n-1
    LR{i} = single(eye(n));
    for k=i+1:n
        LR{i}(k,i) = single(-1);
    end
end

%create difference vectors%
v = cell(n-1);

for k = 1:n-1
    v{k}=single(LR{k}*x);
end


%Create n-1 left L matrices%
LL=cell(n-1);
for i=1:n-1
    M = single(eye(n));
    for k=i+1:n
        M(k,k) = single(1/(v{i}(k)));
    end
    LL{i} = single(M);
end

%Create n-1 U matrices%
N = single((-1)*x);
U=cell(n-1);
for i=1:n-1
    U{i} = single(eye(n));
    for k=i:n-1
        U{i}(k,k+1)= single(N(i));
    end
end

for i=1:n-1
    a = single(LL{i}*single(LR{i}*a));
end

for i=1:n-1
    a = single(U{n-i}*a);
end

