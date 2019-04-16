function val = plot_vph_CN(a,ht,hx)
v = zeros(100,1);
lambda = (a*ht)/(4*hx);
fnc = @(k)((1i)/(k*ht))*log((1-2*(1i)*lambda*sin(k*hx))/(1+2*(1i)*lambda*sin(k*hx)));
for i=1:100
    v(i) = fnc(i);
end
plot(1:100,real(v));
val = 1;