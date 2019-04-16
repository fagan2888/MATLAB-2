function K = ComplexMoment(w,n)

K = zeros(n);

for k = 1:n
    for j = 1:n
        fun = @(t) exp(1i*k*t)*exp(1i*(-j)*t)*w;
        K(k,j) = integral(@(t) fun(t),-pi,pi);
    end
end