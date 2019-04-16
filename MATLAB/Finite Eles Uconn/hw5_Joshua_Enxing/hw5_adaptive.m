%Joshua Enxing
%Hw5 Adaptive Grid Refinement

Nodes = [0 0;1 1;0 1;-1 1];
Triangles = [1 2 3;1 3 4];
T = MakeMesh1(Nodes, Triangles,'Dirichlet');
T.EdgeCFlags(1) = 1;
T.EdgeCFlags(5) = 1;
T.BndyFcn = 'parabolaf';

f=vectorize(inline('exp(-x*y)*(-x^4-2*x^2*y^2-x^2+9*exp(x*y)*(x^4+2*x^2*y^2+x^2+y^4+y^2)*sin(3*x*y)+4*x*y-12*x*y*exp(x*y)*cos(3*x*y)-y^4-y^2)','x','y'));
a=vectorize(inline('x^2+y^2+1','x','y'));
u=vectorize(inline('exp(-x*y)+sin(3*x*y)','x','y'));
ux=vectorize(inline('-y*exp(-x*y)+3*y*cos(3*x*y)','x','y'));
uy=vectorize(inline('-x*exp(-x*y)+3*x*cos(3*x*y)','x','y'));

sol.fnu=u;
sol.fnux=ux;
sol.fnuy=uy;

method = 6;

[U,g,T]=Solve(method,T,a,f,u,[],0.5,50000,2,sol);