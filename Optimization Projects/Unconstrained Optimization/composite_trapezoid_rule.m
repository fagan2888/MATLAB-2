%Joshua Enxing
%Tufts University
%MA226

%Composite Trapezoid Rule
%Takes [a,b] interval of integration,
%function f to integrate (e.g. of the form f = @(x) x^2),
%and number of subpanels N
function int = composite_trapezoid_rule(a,b,N,f)
int = 0;

%Make array of function values at nodes so as to not have to 
%re-compute values each iteration
fn_vals = zeros(N+1,1);
fn_vals(1) = feval(f,a);

%Determine h value
h = (b-a)/N;

%Compute integral contribution from each panel
for i=1:N
    fn_vals(i+1) = feval(f,a+i*h);
    int = int + 0.5*h*(fn_vals(i)+fn_vals(i+1));
end
    