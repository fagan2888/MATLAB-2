%Joshua Enxing
%Code written to approximate the solution to the given problem with 
%lambda = 1 and 5 uniform refinements

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

N=8;

[T,U,g] = TestConv1(T,N,a,f,1,0,u,ux,uy,1);