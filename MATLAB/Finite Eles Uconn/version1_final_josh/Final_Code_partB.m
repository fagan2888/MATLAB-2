%Joshua Enxing
%MA5520
%Final Project Part B

T = RectangleMeshTopLeftD1(1,1,1,1);

%Change the mesh data structure so the nodes and edges are constrained 
%appropriately and all the pointers are correct
T.FNodePtrs = 4;
T.NodePtrs = [-1 -2 -3 1]';
T.CNodePtrs = [1 2 3]';
T.FBndyEdges = [4 5]';
T.EdgeEls(4,2) = -1;
T.EdgeEls(5,2) = -2;
T.EdgeEls(1,2) = 0;

for i=1:6
    T = Refine1(T);
end

%Compute approximate value of h
num_eles = length(T.Nodes(:,1));
h = sqrt(4/num_eles);

t = 0;
dt = 0.01;
k = h;

%We don't know u, but we assign Dirichlet data by making u a function we
%know has the values we want on the Dirichlet part of the boundary. It
%didnt appear to make any difference whether the boundary value at the
%origin was 0 or 1
u = vectorize(inline('1-sign(x)','x','y','t'));
g = getDirichletData(T,u,t);

K = Stiffness1(T,k);
M = Mass(T);
U = zeros(length(M(:,1)),1);
A = Advection(T);
F=Load_New(T,[],k,g,[],t);

for i=1:1/dt
    U = (M+dt*K+dt*A)\(M*U+dt*F);

    figure(1)
    clf
    ShowPWLinFcn1(T,U,g)
    
    t = t + dt;
end

