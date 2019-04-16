function r = residual_source_loc(x,d,A)

r = zeros(5,1);
for i=1:5
    r(i) = (norm(x-A(:,i)).^2-d(i)^2);
end



