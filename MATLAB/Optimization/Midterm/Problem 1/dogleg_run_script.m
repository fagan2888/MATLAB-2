%Joshua Enxing
%Tufts University
%MA150

%Script to run the dogleg method
syms x y;

f = 100*(y-x.^2).^2+(1-x).^2;

trust_region_dogleg(f,[x y],[1.2 1.2]',1,1/4,100,1E-6);
disp(" ");
trust_region_dogleg(f,[x y],[1.2 1.2]',1,1/8,100,1E-6);
disp(" ");
trust_region_dogleg(f,[x y],[1.2 1.2]',1,0,100,1E-6);

disp(" ");

trust_region_dogleg(f,[x y],[-1.2 1]',1,1/4,100,1E-6);
disp(" ");
trust_region_dogleg(f,[x y],[-1.2 1]',1,1/8,100,1E-6);
disp(" ");
trust_region_dogleg(f,[x y],[-1.2 1]',1,0,100,1E-6);