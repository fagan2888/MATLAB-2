x = linspace(-5,5);

y1 = (2 - 1.002.*x)/0.998;
y2 = (-2 - 0.998.*x)/1.002;

figure
plot(x,y1,x,y2);
legend('Line 1','Line 2');
title('6c');

A = [1.002 0.998;0.998 1.002];

kappa = cond(A)

b = [2 -2]';

x = A\b;

e = zeros(1000,2);
r = zeros(1000,2);

for i=1:1000
    deltaA = 0.0001*randn(2,2);
    xi = (A+deltaA)\b;
    e(i,:) = xi' - x';
    r(i,:) = b' - (A*xi)';
end

figure
plot(e(:,1),e(:,2),'.','markersize',20);
title('Forward Error');

figure
plot(r(:,1),r(:,2),'.','markersize',20);
title('Backward Error');
