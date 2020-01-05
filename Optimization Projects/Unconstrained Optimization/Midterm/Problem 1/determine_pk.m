%Joshua Enxing
%Tufts University
%MA150

%Determines pk for the dogleg method
function pk = determine_pk(g,h_mat,pkqn,r)

if norm(pkqn) <= r
        pk = pkqn;
else
    pksd = -(((g')*g)/((g')*h_mat*g))*g;
    pm = pkqn-pksd;
    
    a = (pm')*pm;
    b = 2*(pksd')*pm - 2*(pm')*pm;
    c = (pm')*pm - 2*(pksd')*pm + (pksd')*pksd - r^2;
    
    p = [a b c];
    sols = roots(p);
    t = max(sols);
    
    if 0<=t && t<= 1
        pk = t*pksd;
    else
        pk = pksd + (t-1)*pm;
    end
end