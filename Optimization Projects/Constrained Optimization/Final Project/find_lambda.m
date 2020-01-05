%Finds stepsize using sufficient decrease condition for projected newton
function lambda = find_lambda(alpha,df,beta,f,x,d,u,l)

fxk = f(x(1),x(2));
lambda = 1;

new_x = project(x+lambda*d,u,l);

%check condition, and decrease lambda until met
while f(new_x(1),new_x(2)) - fxk > -alpha*(df')*(x-new_x)
    lambda = beta*lambda;
    new_x = project(x+lambda*d,u,l);
end