%Joshua Enxing%
%University of Connecticut%
%MA5510%

function s = RelativeError(x,f)

[a,b] = leja(x,f);

s = norm((BjorckPereyra(x,f)-BjorckPereyra(a,b))/BjorckPereyra(a,b));