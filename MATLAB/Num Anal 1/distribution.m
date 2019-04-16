function M = distribution(a)
M = zeros(1,ceil(length(a)/2));
p = perms(a);

for i=1:length(p)
    c = 1;
    for j=2:length(a)-1
        b = sort(p(i,1:j));
        b = diff(b);
        d = 1;
        for k=1:length(b)
            if b(k) > 1
                d = d+1;
            end
        end
        c = max(c,d);
    end
    M(c) = M(c)+1;
end