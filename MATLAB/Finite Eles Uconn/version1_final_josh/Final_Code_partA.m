%Joshua Enxing
%MA5520
%Final Project Part A

T = CoarseCircleMeshD1;

k = vectorize(inline('x^2+y^2+1','x','y'));
f = vectorize(inline('6*x*y*(cos(6*x*y*t)-sin(6*x*y*t))+12*t*((3*t*(x^4+2*x^2*y^2+x^2+y^4+y^2)+2*x*y)*sin(6*t*x*y)+(3*t*(x^4+2*x^2*y^2+x^2+y^4+y^2)-2*x*y)*cos(6*t*x*y))','x','y','t'));
u = vectorize(inline('sin(6*x*y*t)+cos(6*x*y*t)','x','y','t'));
ux = vectorize(inline('6*t*y*(cos(6*t*x*y)-sin(6*t*x*y))','x','y','t'));
uy = vectorize(inline('6*t*x*(cos(6*t*x*y)-sin(6*t*x*y))','x','y','t'));

for i=1:4
    T = Refine1(T);
end

t = 0;
dt = 0.01;

K = Stiffness1(T,k);
M = Mass(T);
U = ones(length(M(:,1)),1);

err = 0;

for i=1:(1/dt)

    g = getDirichletData(T,u,t);

    F_next=Load1(T,f,k,g,[],t+dt);

    U = (M+dt*K)\(M*U+dt*F_next);

    figure(1)
    clf
    ShowPWLinFcn1(T,U,g)
    
    err = max(err,L2NormErr1_New(T,u,U,g,t));
    
    t = t + dt;
end

%Compute approximate h by assuming each element is of the form a, a,
%asqrt(2) and use hypotenuse as h
num_eles = length(T.Elements(:,1));
h = sqrt((4*pi)/num_eles);

disp('LinfL^2 err       dt           h');
disp([err               dt           h]);


    
    