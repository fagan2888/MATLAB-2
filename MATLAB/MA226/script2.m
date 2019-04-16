%Joshua Enxing
%Tufts University
%MA226 

%Script for HW 7 3(c)

g2 = @(x)sqrt(x+2);
g3 = @(x)1+2./x;
g4 = @(x)(x.^2+2)/(2.*x-1);

e2 = zeros(100,1);
e3 = zeros(100,1);
e4 = zeros(100,1);

g20 = 1;
g30 = 1;
g40 = 1;

g2stop = 100;
g3stop = 100;
g4stop = 100;

for i=1:100
    e2(i) = abs(2-g20);
    e3(i) = abs(2-g30);
    e4(i) = abs(2-g40);
    
    if e2(i) > 10^(-6)
        g20 = g2(g20);
    else
        g2stop = min(g2stop,i);
    end
    
    if e3(i) > 10^(-6)
        g30 = g3(g30);
    else
        g3stop = min(g3stop,i);
    end
    
    
    if e4(i) > 10^(-6)
        g40 = g4(g40);
    else
        g4stop = min(g4stop,i);
    end
    
    if e2(i) <= 10^(-6) && e3(i) <= 10^(-6) && e4(i) <= 10^(-6)
        break;
    end
    
end

e2 = e2(1:g2stop);
e3 = e3(1:g3stop);
e4 = e4(1:g4stop);

l2x = log(e2);
l2x = l2x(1:length(l2x)-1);
l2y = log(e2);
l2y = l2y(2:length(l2y));

l3x = log(e3);
l3x = l3x(1:length(l3x)-1);
l3y = log(e3);
l3y = l3y(2:length(l3y));

l4x = log(e4);
l4x = l4x(1:length(l4x)-1);
l4y = log(e4);
l4y = l4y(2:length(l4y));

plot(l2x,l2y,l3x,l3y,l4x,l4y)
legend('Log error g_2(x)','Log error g_3(x)','Log error g_4(x)')
xlabel('log(e_k)');
ylabel('log(e_{k+1})');
title('Log(e_k) versus Log(e_{k+1})');
    
