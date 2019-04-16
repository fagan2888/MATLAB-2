function x = make_input_vector(k,n)

x = zeros(n,1);

for i=1:n
    x(i) = sin((k*pi*i)/n+1);
end

