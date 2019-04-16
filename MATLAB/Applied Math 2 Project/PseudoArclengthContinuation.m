%Joshua Enxing
%University of Connecticut
%MA5411
%12/7/2015

function [Iu, L] = PseudoArclengthContinuation(ds,p,N)

%Fill vector for numerical integration
t = zeros(N-1,1);
for p=1:N-1
    t(p) = (1/N)*p;
end

%Step size is 1/(#subintervals in mesh)
h = 1/N;
u = zeros(p,N-1);
tanv = zeros(p,N);
for i=1:N
    tanv(1,i) = 1;
end
tanv(1,:) = tanv(1,:)/norm(tanv(1,:));
dul = zeros(N,1);
v = zeros(N,1);
v(N) = 1;
    
%lambda_min is 0 here
l = 0;
s = 0

M = zeros(N);
R = zeros(N,1);

%F is F(u,lambda) where F(k) = F(u_k, lambda)	
F = zeros(N-1,1);

%Fu is F_u
Fu = zeros(N-1);

%Iu(k) is the integral of each u_k
Iu = zeros(p,1);
L = zeros(p,1);

%Fill super- and sub-diagonals of F_u
for i=1:N-2
    Fu(i,i+1) = 1/h^2;
    Fu(i+1,i) = 1/h^2;
end

%G = F_lambda
G = zeros(N-1,1);

for k = 1:p
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
    
    %Fill G = F_lambda
    for r=1:N-1
        G(r) = exp(u(k,r));
    end
    
    if k==1
        M(1:N-1,1:N-1) = Fu;
        M(1:N-1,N) = G;
        M(N,1:N-1) = tanv(k,1:N-1);
        M(N,N) = dul(N);
        if k == 1
            d = dot(u(k,:),tanv(k,1:N-1)') - ds;
        else
            d = dot(u(k,:)-u(k-1,:),tanv(k-1,1:N-1)')+(L(k)-L(k-1))*tanv(k-1,N) - ds;
        end
        R(1:N-1,1) = F';
        R(N,1) = d;
        R = (-1)*R;
        dul = M\R;
        
        u(k,:) = u(k,:) + dul(1:N-1)';
        l = l + dul(N);
        L(k) = l;
 
    
        %Update F
        for i=2:N-2
            F(i) = (u(k,i+1) - 2*u(k,i) + u(k,i-1))/h^2 + l*exp(u(k,i));
        end
    
        F(N-1) = l*exp(u(k,N-1)) - (2*u(k,N-1))/h^2 + u(k,N-2)/h^2;
        F(1) = (u(k,2) - 2*u(k,1))/h^2 + l*exp(u(k,1));
        
        for r=1:N-1
            G(r) = exp(u(k,r));
        end
    end
    
    %Newton's Method until converges
    while (norm(abs(F)) > sqrt(eps)) 
        M(1:N-1,1:N-1) = Fu;
        M(1:N-1,N) = G;
        M(N,1:N-1) = tanv(k,1:N-1)';
        M(N,N) = dul(N);
        if k == 1
            d = dot(u(k,:),tanv(k,1:N-1)') - ds;
        else
            d = dot(u(k,:)-u(k-1,:),tanv(k-1,1:N-1)')+(L(k)-L(k-1))*tanv(k-1,N) - ds;
        end
        R(1:N-1,1) = F';
        R(N,1) = d;
        R = (-1)*R;
        dul = M\R;
        
        u(k,:) = u(k,:) + dul(1:N-1)';
        l = l + dul(N);
        L(k) = l;
        
        %Update F
        for i=2:N-2
            F(i) = (u(k,i+1) - 2*u(k,i) + u(k,i-1))/h^2 + l*exp(u(k,i));
        end
    
        F(N-1) = l*exp(u(k,N-1)) - (2*u(k,N-1))/h^2 + u(k,N-2)/h^2;
        F(1) = (u(k,2) - 2*u(k,1))/h^2 + l*exp(u(k,1));
        
        for r=1:N-1
            G(r) = exp(u(k,r));
        end
    end
    
    %Compute next direction vector and normalize
    tanv(k+1,:) = M\v;
    tanv(k+1,:) = tanv(k+1,:)/norm(tanv(k+1,:),2);
    
    %Numerically integrate u
    Iu(k) = trapz(t, u(k,:));
end

%Graph numerical approximation of bifurcation diagram
%plot(L,Iu);
