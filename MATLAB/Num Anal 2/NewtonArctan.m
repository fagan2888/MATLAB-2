function a = NewtonArctan(x)
a = x - atan(x)*(1+x^2);
b = x;
while abs(b-a) > 10^(-8) 
    b = a;
    a = a - atan(a)*(1+a^2);
    disp(b);
    disp(a);
end


    