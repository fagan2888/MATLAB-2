function [D,x] = make_diff_matrix(N)


x = make_collocation_points(N);
D = zeros(length(x));

for i=2:N-1
    if i > 1 && i < N
        D(i,i) = -x(i)/(2*(1-(x(i))^2));
    end
end

D(1,1) = (2*(N-1)^2+1)/6;
D(N,N) = -1*D(1,1);
for i=1:N
    for j=1:N
        if i==1 || i==N
            ci = 2;
        else
            ci = 1;
        end
        if j==1 || j==N
            cj = 2;
        else
            cj = 1;
        end
        if i ~= j
            D(i,j) = ((ci/cj)*((-1)^(i+j)))/(x(i)-x(j));
        end
    end
end