%Evaluates the boundary condition functions at all values of t and
%store these values in vectors
function [l_bcs,r_bcs] = make_bc_vecs(left_bcs,right_bcs,t_nodes,h_t)

l_bcs = zeros(t_nodes,1);
r_bcs = zeros(t_nodes,1);

for i=1:t_nodes
    l_bcs(i) = left_bcs((i-1)*h_t);
    r_bcs(i) = right_bcs((i-1)*h_t);
end
