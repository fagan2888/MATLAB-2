%Joshua Enxing
%University of Connecticut
%MA5511

function z = Maehly(f)
%Find f'(x)
k = diff(f);

%Create array with coefficients of f, also find degree of f
A = zeros([1,1]);
j = 1;
if abs(subs(f,0)) > 0
    A(1) = subs(f,0);
    j = j+1;
end
r = diff(f);
i = 2;
while r ~= 0
    t = subs(r,0)/factorial(i-1);
    if (t ~= 0)
        A(j) = t;
        j = j+1;
    end
    i = i+1;
    r = diff(r);
end
n = i-2;
c = A;
    

%Create array to hold bound on roots
X = zeros([1,5]);

%1st bound on roots
A = zeros([1,length(c)-1]);
A(1) = abs(c(1)/length(c));
for i = 2:length(c)-1
    A(i) = 1 + abs(c(i)/c(length(c)));
end
a = max(A);
X(1) = a;

%2nd bound on roots
B = zeros([1,2]);
b = 0;
for i = 1:length(c)-1
    b = b + abs(c(i)/c(length(c)));
end
B(1) = 1;
B(2) = b;
b = max(B);
X(2) = b;

%3rd bound on roots
D = zeros([1,length(c)-1]);
D(1) = abs(c(1)/c(2));
for i = 2:length(c)-1
    D(i) = 2*abs(c(i)/c(i+1));
end
d = max(D);
X(3) = d;

%4th bound on roots
e = 0;
for i = 1:length(c)-1
    e = e + abs(c(length(c)-i)/c(length(c)-i+1));
end
X(4) = e;

%5th bound on roots
G = zeros([1,length(c)-1]);
G(1) = abs(c(length(c)-1)/c(length(c)));
for i = 2:length(c)-1
    G(i) = nthroot(abs(c(length(c)-i)/c(length(c))),i);
end
g = 2*max(G);
X(5) = g;

%Lowest of all computed bounds on roots
x = min(X);

%Create array to hold computed roots
z = zeros([1,n]);

%Make l arb large so we will pass first check of while loop
h = x;
l = x+2;

%One iteration of j for each degree of f.
for j = 1:n
    %Double-step
    m = 2;
    while(abs(l-x) > 10^(-8))
        s = 0;
        l = h;
        h = x;
        
        %Make sum
        for i = 1:j-1
            s = s + 1/(h-z(i));
        end
        
        %Compute next iterate
        x = h - (m*subs(f,h)/(subs(k,h)-s*subs(f,h)));
        
        %Do Newton's Method (single-step) with previous iterate if overshoot
        if (h < x && m == 2)
            m = 1;
            x = h;
        end
    end
    
    %Found root
    z(j) = x;
    
    %Make next guess less than root, since we find largest root first
    x = x - 10^(-4);
end
