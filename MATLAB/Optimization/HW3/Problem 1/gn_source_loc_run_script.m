 rng(317,'v4');
 A = randn(2,5); 
 x = randn(2,1); 
 d = sqrt(sum((A-x*ones(1,5)).^2) )+0.05*randn(1,5);
 
 [rs,z] = gauss_newton_source_loc([2 2]',d,A,1,0.9,0.5,100,1E-6);
 
 figure
 plot(rs);
 xlabel('Iteration');
 ylabel('Norm of residual');
 title('x_0 = [2 2]^T');
 
 [rs,z] = gauss_newton_source_loc([1 1]',d,A,1,0.9,0.5,100,1E-6);
 
 
 figure
 plot(rs);
 xlabel('Iteration');
 ylabel('Norm of residual');
 title('x_0 = [1 1]^T');
 
 [rs,z] = gauss_newton_source_loc([0.5 0]',d,A,1,0.9,0.5,100,1E-6);
 
 
 figure
 plot(rs);
 xlabel('Iteration');
 ylabel('Norm of residual');
 title('x_0 = [0.5 0]^T');