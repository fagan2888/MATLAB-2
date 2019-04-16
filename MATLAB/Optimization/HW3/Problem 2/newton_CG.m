function x = newton_CG(x0,nk,alpha0,c,rho,max_iters,tol,option)

x = x0;

disp("Iteration || Norm Gradient || nk");
for k=1:max_iters
    fx = extended_rosenbrock(x);
    gradf = grad_f(x);
    hessf = hess_f(x);
    
    disp(k-1 + " || " + norm(gradf) + " || " + nk);
    if norm(gradf) < tol
        break;
    end
    
    if option==1
        nk = min(0.5,sqrt(norm(gradf)));
    end
    if option==2
        nk = min(0.5,norm(gradf));
    end
        
    pk = CG_inexact_Newton(hessf,gradf,nk,max_iters);
    
    alpha = alpha0;
    while extended_rosenbrock(x+alpha*pk) > fx + c*alpha*(gradf')*pk
       alpha = rho*alpha;
    end
    %alpha = find_alpha(x,pk,0,0.9,1,fx,gradf);
    
    x = x + alpha*pk;
end