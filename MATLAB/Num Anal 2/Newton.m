%Joshua Enxing
%University of Connecticut
%MA5511

function a = Newton(f,x)
%Find f'(x)
g = diff(f);

%Find first iterate
a = x - subs(f,x)/subs(g,x);
b = x;

%Keep track of number of iterations
i = 1;

%While difference between iterates is greater than 10^(-8)
%and while number of iterates is less than 100
%then we keep iterating

%Make header for output of iterates and the ratios
S = {'Iterate' 'Ratio'};
disp(S);
while abs(b-a) > 10^(-8) && i < 100
    %c is (k-1)st iterate
    c = b;
    %b is kth iterate
    b = a;
    %a is (k+1)st iterate
    a = b - subs(f,b)/subs(g,b);
    
    %Increment iterate counter
    i = i+1; 
    
    d = abs(a-b);
    e = abs(b-c)^2;
    S = {num2str(b) num2str(d/e)};
    disp(S);
end


    