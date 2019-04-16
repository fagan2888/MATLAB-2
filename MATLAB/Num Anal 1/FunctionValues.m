function [x,f] = FunctionValues

x = zeros(14,1);
f = zeros(14,1);
x(1) = -10;
x(2) = -8;
x(3) = -7;
x(4) = -5;
x(5) = -3;
x(6) = -2;
x(7) = -1/2;
x(8) = 1;
x(9) = 2;
x(10) = 3;
x(11) = 5;
x(12) = 7;
x(13) = 8;
x(14) = 10;
for i=1:14
    f(i) = x(i)^13;
end

