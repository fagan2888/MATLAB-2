%Joshua Enxing
%University of Connecticut
%MA5411
%12/7/2015

function [Iu, L] = ParameterContinuation(dl,s,N)

%Fill vector for numerical integration
t = zeros(N-1,1);
for p=1:N-1
    t(p) = (1/N)*p;
end

%Step size is 1/(#subintervals in mesh)
h = 1/N;
u = zeros(s,N-1);

%lambda_min is 0 here
l = 0;
du = zeros(N-1,1);

%F is F(u,lambda) where F(k) = F(u_k, lambda)	
F = zeros(N-1,1);

%Fu is F_u
Fu = zeros(N-1);

%Iu(k) is the integral of each u_k
Iu = zeros(s,1);
L = zeros(s,1);

%Fill super- and sub-diagonals of F_u
for i=1:N-2
    Fu(i,i+1) = 1/h^2;
    Fu(i+1,i) = 1/h^2;
end

%G = F_lambda
G = zeros(N-1,1);

for k = 1:s
    
    if k > 1
	%Start with initial guess for u_k
        u(k,:) = u(k-1,:) + dl*(du');
    end

    %Increment lambda
    l = l + dl;
    L(k) = l;
    
    for i=1:N-1
	%Fill diagonal entries of F_u
        Fu(i,i) = l*exp(u(k,i))-2/h^2;
    end
         
    %Fill F
    for i=2:N-2
        F(i) = (u(k,i+1) - 2*u(k,i) + u(k,i-1))/h^2 + l*exp(u(k,i));
    end
    
    F(N-1) = l*exp(u(k,N-1)) - (2*u(k,N-1))/h^2 + u(k,N-2)/h^2;
    F(1) = (u(k,2) - 2*u(k,1))/h^2 + l*exp(u(k,1));
    
    %Newton's Method until converges
    while (norm(abs(F)) > sqrt(eps))
        F = F*(-1);
        x = Fu\F;
        F = F*(-1);
        u(k,:) = u(k,:) + x';
        
	%Update F
        for i=2:N-2
            F(i) = (u(k,i+1) - 2*u(k,i) + u(k,i-1))/h^2 + l*exp(u(k,i));
        end
    
        F(N-1) = l*exp(u(k,N-1)) - (2*u(k,N-1))/h^2 + u(k,N-2)/h^2;
        F(1) = (u(k,2) - 2*u(k,1))/h^2 + l*exp(u(k,1));
    end
    
    %Numerically integrate u
    Iu(k) = trapz(t, u(k,:));
    
    %Fill G = F_lambda
    for r=1:N-1
        G(r) = exp(u(k,r));
    end
    G = G*(-1);
    
    %Find tangent vector so can use for guess in next iterate
    du = Fu\G;
    G = G*(-1);
end

%Graph numerical approximation of bifurcation diagram
plot(L,Iu);




