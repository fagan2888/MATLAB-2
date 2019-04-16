function err_vect = produce_error(eps_grid)

a = rand();
b = rand();

n = length(eps_grid);
err_vect = zeros(1,n);

for i = 1:n
    err_vect(i) = err_for_eps(eps_grid(i),a,b);
    disp(i);
end

loglog(eps_grid,err_vect);