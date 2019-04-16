g1 = @(x)x.^2-2;
g2 = @(x)sqrt(x+2);
g3 = @(x)1+2./x;
g4 = @(x)(x.^2+2)/(2.*x-1);

x = ones(4,1);
y = ones(4,1);
z = ones(4,1);
w = ones(4,1);

for i=2:4
    x(i) = g1(x(i-1));
    y(i) = g2(y(i-1));
    z(i) = g3(z(i-1));
    w(i) = g4(w(i-1));
end

a = [0 1 2 3]';

figure
plot(a,x,a,y,a,z,a,w);
legend('g_1(x) = x^2-2','g_2(x) = sqrt(x+2)','g_3(x) = 1+2/x','g_4(x) = (x^2+2)/(2x-1)');
xlabel('Iterate Number');
ylabel('Iterate Value');
title('x_0 = 1');

x = 2.5*ones(4,1);
y = 2.5*ones(4,1);
z = 2.5*ones(4,1);
w = 2.5*ones(4,1);

for i=2:4
    x(i) = g1(x(i-1));
    y(i) = g2(y(i-1));
    z(i) = g3(z(i-1));
    w(i) = g4(w(i-1));
end

figure
plot(a,x,a,y,a,z,a,w);
legend('g_1(x) = x^2-2','g_2(x) = sqrt(x+2)','g_3(x) = 1+2/x','g_4(x) = (x^2+2)/(2x-1)');
xlabel('Iterate Number');
ylabel('Iterate Value');
title('x_0 = 2.5');
