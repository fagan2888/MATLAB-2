function x = lbfgs_rosen(x,m,alpha0,rho,c,max_iters,tol)

S = zeros(length(x),m);
Y = zeros(length(x),m);
gfxk = grad_f(x);
p = zeros(m,1);

disp("Iteration || norm(grad)");
for k=1:max_iters
    fxk = extended_rosenbrock(x);
    
    format long;
    disp(k-1 + " || " + norm(gfxk));
    if norm(gfxk) < tol
        break;
    end
    
    if k == 1
        pk = -gfxk;
    else
        pk = compute_pk(p,gfxk,S,Y,k-1,m);
    end
    
    alpha = alpha0;
    fxpk = extended_rosenbrock(x+alpha*pk);
    while fxpk > fxk + c*alpha*(gfxk')*pk
       alpha = rho*alpha;
       fxpk = extended_rosenbrock(x+alpha*pk);
    end
    
    %alpha = find_alpha(x,pk,0,0.8,0.9,fxk,gfxk);
    
    xlast = x;
    x = x + alpha*pk;
    gfxk_last = gfxk;
    gfxk = grad_f(x);
    
    if k > m
        S(:,k-m) = zeros(length(x),1);
        Y(:,k-m) = zeros(length(x),1);
        p(k-m) = 0;
    end
    S(:,k) = x - xlast;
    Y(:,k) = gfxk - gfxk_last;
    p(k) = 1/((Y(:,k)')*S(:,k));
%     disp(S);
%     disp(Y);
%     disp(p);
    
end
    
    
    
    