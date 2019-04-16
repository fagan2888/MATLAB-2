function f = Untitled2(n)

LR = cell(n-1);

for i=1:n-1
    LR{i} = eye(n);
    for k=i+1:n
        LR{i}(k,i) = -1;
    end
end

for l=1:n-1
    disp(LR{l});
end
