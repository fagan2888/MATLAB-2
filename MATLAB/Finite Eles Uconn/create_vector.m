%creates vector v for an n-element partion of [a,b]
function v = create_vector(a,b,m)

h = (b-a)/(m+1);

v = zeros(m+2,1);

for i=1:m+2
    v(i) = (i-1)*h + a;
end

