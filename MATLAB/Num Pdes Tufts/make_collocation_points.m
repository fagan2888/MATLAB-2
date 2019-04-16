function x = make_collocation_points(N)

x = zeros(N,1);

for i=1:N
    x(i) = cos((pi*(i-1))/(N-1));
end
