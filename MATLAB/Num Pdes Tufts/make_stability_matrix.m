function [R,dt] = make_stability_matrix(D2,dt)
disp(dt);
N = length(D2(:,1));

R = eye(N) + dt*D2 + 0.5*dt^2*(D2^2) + (dt^3/6)*(D2^3) + (dt^4/24)*(D2^4);

% d = eigs(R);
% disp(d);
% 
% m = norm(d(1));
% for i=2:6
%     d(i) = norm(d(i));
%     m = max(m,d(i));
% end
% 
% 
% if m > 1
%     [R,dt] = make_stability_matrix(D2,dt/2);
% end