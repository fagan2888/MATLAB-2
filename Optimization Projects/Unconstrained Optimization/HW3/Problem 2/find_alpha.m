function alpha = find_alpha(x,pk,a1,a2,amax,fxk,gradf)

i = 1;
for j=1:100
    fa2 = extended_rosenbrock(x+a2*pk);
    fa1 = extended_rosenbrock(x+a1*pk);
    if (fa2 > fxk + (1/3)*a2*(pk')*gradf) || (fa2 >= fa1 && i > 1)
        alpha = zoom(x,pk,a1,a2,fxk,gradf);
        break;
    end
    gfa2 = grad_f(x+a2*pk);
    if abs((pk')*gfa2) <= -(2/3)*(pk')*gradf
        alpha = a2;
        break;
    end
    if (pk')*gfa2 >= 0
        alpha = zoom(x,pk,a2,a1,fxk,gradf);
        break;
    end
    a1 = a2;
    a2 = 0.5*(amax+a2);
    i = i+1;
end
    

    

