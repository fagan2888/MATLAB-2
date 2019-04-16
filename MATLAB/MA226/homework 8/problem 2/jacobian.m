 function J = jacobian(w0,w1,x0,x1)

J = zeros(4);
J(1,:) = [1 1 0 0];
J(2,:) = [x0 x1 w0 w1];
J(3,:) = [x0^2 x1^2 2*w0*x0 2*w1*x1];
J(4,:) = [x0^3 x1^3 3*w0*x0^2 3*w1*x1^2];