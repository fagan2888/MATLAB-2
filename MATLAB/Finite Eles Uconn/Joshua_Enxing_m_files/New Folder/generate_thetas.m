%Joshua Enxing
%University of Connecticut
%MA5520
%Due 2/19/16

%generates the piece-wise linear finite elements for this problem (note
%that the elements are defined past their support for convenience, but they
%are only integrated on their support
function theta = generate_thetas(z, n)
syms x;
theta = sym('theta%d%d',[n+2 1]);

theta(1) = sym((x-z(2))/(0-z(2)));
theta(n+2) = sym((x-z(n+1))/(1-z(n+1)));

for i=2:n+1
    theta(i) = sym(heaviside(z(i)-x)*((x-z(i-1))/(z(i)-z(i-1))) + heaviside(x-z(i))*((x-z(i+1))/(z(i)-z(i+1))));
end
