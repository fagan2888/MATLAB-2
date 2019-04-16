function f = ThirteenthPower(x)

n = length(x);
f = zeros(n,1);
for i=1:n
    f(i) = x(i)^13;
end

