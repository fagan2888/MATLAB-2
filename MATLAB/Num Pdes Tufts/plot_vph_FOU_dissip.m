function val = plot_vph_FOU_dissip(a,ht,hx,l)
v = zeros(100,1);
lambda = (a*ht)/hx;
fnc = @(k)abs(1-lambda+lambda*exp(-1*(1i)*k*hx))^l;
for i=1:100
    v(i) = fnc(i);
end
plot(1:100,real(v));
val = 1;