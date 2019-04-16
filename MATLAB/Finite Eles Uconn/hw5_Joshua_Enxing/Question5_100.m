%Joshua Enxing
%Code written to approximate the solution to the given problem with 
%lambda = 100 and 5 uniform refinements

Nodes = [0 0;1 1;0 1;-1 1];
Triangles = [1 2 3;1 3 4];
T = MakeMesh1(Nodes, Triangles,'Dirichlet');
T.EdgeCFlags(1) = 1;
T.EdgeCFlags(5) = 1;
T.BndyFcn = 'parabolaf';

f=vectorize(inline('-(2*x*(-100*y*exp(-100*x*y)+3*y*cos(3*x*y))+(x^2+y^2+1)*(10000*y^2*exp(-100*x*y)-9*y^2*sin(3*x*y))+2*y*(-100*x*exp(-100*x*y)+3*x*cos(3*x*y))+(x^2+y^2+1)*(10000*x^2*exp(-100*x*y)-9*x^2*sin(3*x*y)))','x','y'));
a=vectorize(inline('x^2+y^2+1','x','y'));
u=vectorize(inline('exp(-100*x*y)+sin(3*x*y)','x','y'));
ux=vectorize(inline('-100*y*exp(-100*x*y)+3*y*cos(3*x*y)','x','y'));
uy=vectorize(inline('-100*x*exp(-100*x*y)+3*x*cos(3*x*y)','x','y'));

N=5;

[T,U,g] = TestConv1(T,N,a,f,1,0,u,ux,uy,1);