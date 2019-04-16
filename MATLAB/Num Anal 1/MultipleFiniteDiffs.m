%Joshua Enxing
%University of Connecticut
%MA5511

function MultipleFiniteDiffs

FE = zeros(1,16);
BE = zeros(1,16);
MID = zeros(1,16);
CN = zeros(1,16);
FE(1:16) = 10;
BE(1:16) = 10;
MID(1:16) = 10;
CN(1:16) = 10;

for k = 2:9
    dt = 1/2^k;
    for j = 1:2^k
        FE(k-1) = FE(k-1) + dt*(-1)*FE(k-1);
        BE(k-1) = BE(k-1)/(1-(-1)*dt);
        MID(k-1) = MID(k-1)*(1+dt*(1+(dt/2)))*(-1);
        CN(k-1) = ((1+(dt/2)*(-1))/(1-(dt/2)*(-1)))*CN(k-1);
        
        FE(2*k-2) = FE(2*k-2) + dt*(-50)*FE(2*k-2);
        BE(2*k-2) = BE(2*k-2)/(1-(-50)*dt);
        MID(2*k-2) = MID(2*k-2)*(1+dt*(1+(dt/2)))*(-50);
        CN(2*k-2) = ((1+(dt/2)*(-50))/(1-(dt/2)*(-50)))*CN(2*k-2);
    end
end

X = zeros(1,8);
for k = 2:9
    X(k-1) = log(1/2^k);
    FE(k-1) = log(abs((1/exp(1))-FE(k-1)));
    BE(k-1) = log(abs((1/exp(1))-BE(k-1)));
    MID(k-1) = log(abs((1/exp(1))-MID(k-1)));
    CN(k-1) = log(abs((1/exp(1))-CN(k-1)));
    
    FE(2*k-2) = log(abs((1/exp(50))-FE(2*k-2)));
    BE(2*k-2) = log(abs((1/exp(50))-BE(2*k-2)));
    MID(2*k-2) = log(abs((1/exp(50))-MID(2*k-2)));
    CN(2*k-2) = log(abs((1/exp(50))-CN(2*k-2)));
end

plot(X,FE(1:8),X,BE(1:8),X,MID(1:8),X,CN(1:8));
%plot(X,FE(9:16),X,BE(9:16),X,MID(9:16),X,CN(9:16));




