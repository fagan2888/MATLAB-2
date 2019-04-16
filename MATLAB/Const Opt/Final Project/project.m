%Projects x onto boundary of constraints with upper and lower bounds u and
%l
function x = project(x,u,l)

n = length(x);

for i=1:n
    if x(i) > u(i)
        x(i) = u(i);
    end
    if x(i) < l(i)
        x(i) = l(i);
    end
end