%  function [ model ] = prob_gen
%
%     Purpose:
%
%     Generates matrices and discretization for the diffusion eqn,
%
%                   f_t - theta * f_xx = 0
%
%              and periodic boundary conditions
%
%                   f (0) = f (1),
%                   f'(0) = f'(1).
%
%              with theta(1) the diffusion coefficient for 0 <= x < 1/2,
%                   theta(2) the diffusion coefficient for 1/2 <= x < 1.
%
%  Version Jan 19, 2015
%  Caleb Magruder
%
function [ model ] = prob_gen

    nx = 100; %%%MUST BE EVEN!
    nt = 40;
    tf = 2;
    
    y0 = [ 1*ones(nx/2,1); 0*ones(nx/2,1) ];
    
    dx = 1 / nx;
    dt = tf / nt;
    
    t = linspace(0,tf,nt+1);
    
    %%%Stiffness Matrix
    e1_sub = -ones(nx/2,1);
    e1_diag = [2*ones((nx-2)/2,1); 1];
    A1 = 1/dx*spdiags([e1_sub, e1_diag, e1_sub],[-1:1],nx/2,nx/2);
    A1(nx,nx) = 1/dx;    A1(nx,1) = -1/dx;    A1(1,nx) = -1/dx;
    
    e2_sub = -ones(nx/2 + 1,1);
    e2_diag = [1; 2*ones((nx-2)/2,1); 1];
    A2 = sparse(nx,nx);
    A2(nx/2:nx, nx/2:nx) = ...
        1/dx*spdiags([e2_sub, e2_diag, e2_sub],[-1:1],nx/2+1,nx/2+1);
    
    A = @(k) k(1)*A1 + k(2)*A2;
    
    %%%Mass Matrix
    e = ones(nx, 1);
    M       = (dx/6)*spdiags([e 4*e e], -1:1, nx, nx);
    M(1,nx) = (dx/6);
    M(nx,1) = (dx/6);
    
    p = 2;
    
    model = struct('M',M,'A',A,'A1',A1,'A2',A2,...
        't',t,'dx',dx,'dt',dt,'nx',nx,'nt',nt,'tf',tf,'y0',y0,'p',p);

end