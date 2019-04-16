function val = plot_vph_FOU(a,ht,hx)
v = zeros(100,1);
lambda = (a*ht)/hx;
fnc = @(k)((1i)/(k*ht))*log((1-lambda+lambda*exp(-1*(1i)*k*hx))/abs(1-lambda+lambda*exp(-1*(1i)*k*hx)));
for i=1:100
    v(i) = fnc(i);
end
plot(1:100,real(v));
val = 1;