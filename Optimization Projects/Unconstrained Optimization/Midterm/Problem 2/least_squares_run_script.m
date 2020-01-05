[x,y] = generate_points;
[a,b,c] = quadratic_least_squares(1,1,1,1,0.9,0.5,100,1E-6);

disp("a = " + a + ", b = " + b + ", c = " + c);

z = a.*x.^2 + b.*x + c;
plot(x,y,x,z);
legend('Data Points','Computed Polynomial');
