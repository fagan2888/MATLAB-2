function pk = compute_pk(p,gfxk,S,Y,k,m)

q = gfxk;

alphas = zeros(m,1);
for i=1:m
    j = k-i;
    if j > 0 
        pj = p(j);
        alphas(i) = pj*(S(:,j)')*q;
        q = q - alphas(i)*Y(:,j);
    end
 end

gamma = ((S(:,k)')*Y(:,k))/((Y(:,k)')*Y(:,k));
r = gamma*q;

for i=0:m-1
    j = k-m+i;
    if j > 0
        pj = p(j);
        beta = pj*(Y(:,j)')*r;
        
        r = r + S(:,j)*(alphas(k-j)-beta);
    end
end

pk = -r;
