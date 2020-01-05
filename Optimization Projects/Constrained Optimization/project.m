%Projection for projected gradiend on [-1,1]^2
function z = project(z)

%If first component outside domain
if z(1) > 1
    z(1) = 1;
elseif z(1) < -1
    z(1) = -1;
end

%If second component outside domain
if z(2) > 1
    z(2) = 1;
elseif z(2) < -1
    z(2) = -1;
end