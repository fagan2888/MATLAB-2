%Joshua Enxing
%University of Connecticut
%MA5511

function b = Broyden(f,x)
%Set h of reasonable size
h = sqrt(eps);

%Set a arbitrary so that abs(b-a) > 10^(-8)
a = x+10;
b = x;

%Iterate counter 
j = 0;

%Make B_0
B = (subs(f,a+h)-subs(f,a))/h;

%Make header for output of iterates and the ratios
S = {'Iterate' 'Ratio'};
disp(S);
while abs(b-a) >= 10^(-8) && j < 100
    %Set d_k
    d = subs(f,b)/B;
    i = 0;
    
    %Find i such that |f(x_k-2^(-i)| <= |f(x_k)|
    while abs(subs(f,b-(2^(-i))*d)) > abs(subs(f,b))
        i = i+1;
    end
    
    %Set lambda_k
    l = 2^(-i);
    
    %c is (k-1)st iterate
    c = a;
    
    %a is kth iterate
    a = b;
    
    %b is (k+1)st iterate
    b = a - l*d;
    j = j+1;
    
    %Compute ratios
    d = abs(b-a);
    e = abs(a-c)^2;
    S = {num2str(a) num2str((d/e))};
    
    %Once we have the second iterate we will display the ratios
    if j >= 2
        disp(S);
    end
    
    %Determine B_(k+1)
    if i <= 1
        p = b - a;
        q = subs(f,b) - subs(f,a);
        B = B + (1/p^2)*(q-B*p)*p;
    else
        B = (subs(f,b+h)-subs(f,b))/h;
    end
end
