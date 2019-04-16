function [left_bcs,right_bcs] = make_bc_vec(left_bcs,right_bcs,t_nodes,h_t)

bcs = linspace(0,1,t_nodes);

for i=1:t_nodes
    bcs(i) = bdy_cond((i-1)*h_t);
end
