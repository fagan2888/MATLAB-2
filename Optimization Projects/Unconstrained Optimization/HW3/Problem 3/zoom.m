function alpha = zoom(x,pk,a1,a2,fxk,gradf)

boo = true;
while boo
    ai = 0.5*(a1+a2);
    fai = extended_rosenbrock(x+ai*pk);
    fa1 = extended_rosenbrock(x+a1*pk);
    if (fai > fxk + (1/3)*ai*(pk')*gradf) || (fai >= fa1)
        a2 = ai;
    else
        gfai = grad_f(x+ai*pk);
        if (pk')*gfai <= (2/3)*(pk')*gradf
        	alpha = ai;
            boo = false;
        end
        if ((pk')*gfai)*(a2-a1) >= 0
            a2 = a1;
        end
        a1 = ai;
    end
end