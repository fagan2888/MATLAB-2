%Joshua Enxing
%University of Connecticut
%MA5520
%Due 2/19/16

%calculate the l2 and h1 errors for our FEM solutions for the number of 
%elements in each entry of the vector num_eles, and q is our vector of
%h values corresponding to each mesh
function [q,l2,h1] = calculate_error_vectors(num_eles,u,alpha,beta,a,b,f,m)

syms x;
q = zeros(length(num_eles),1);
l2 = zeros(length(num_eles),1);
h1 = zeros(length(num_eles),1);

for i=1:length(num_eles)
    %make mesh corresponding to the number of elements prescribed in the
    %ith entry of the vector num_eles
    mesh = linspace(0,1,num_eles(i)+1);
    
    %run our FEM for this mesh
    [u_h,theta,h] = pw_linear_FEM_mixed_BCs(mesh,alpha,beta,a,b,f,m);
    
    %put this h value into the q vector
    q(i) = h(1);
    
    %calculate the l2 and h1 errors for FEM solution with this mesh
    [h1(i),l2(i)] = h1_error(u,alpha,beta,u_h,theta,h,mesh,m);
end