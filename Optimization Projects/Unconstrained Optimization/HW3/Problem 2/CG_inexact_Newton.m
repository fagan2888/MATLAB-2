function pk = CG_inexact_Newton(hess,grad,nk,max_iters)

z = 0;
r = grad;
d = -r;

for j=1:max_iters
    if (d')*hess*d < 0
        if j==0
            pk = -grad;
        else
            pk = z;
        end
        break;
    end
    
    alpha = ((r')*r)/((d')*hess*d);
    z = z + alpha*d;
    
    r_last = r;
    r = r + alpha*hess*d;
    
    if norm(r) <= nk*norm(grad)
        pk = z;
        break;
    end
    
    beta = ((r')*r)/((r_last')*r_last);
    d = -r+beta*d;
end