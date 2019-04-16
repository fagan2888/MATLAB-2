function val = plot_L1_solution(a,A0,k,x0,xn,num_nodes,tn,num_t_nodes)

%L1
%{
dt = 1/100;
x_nodes = linspace(0,2,100);
B = @(z)0;
A = @(q)q;
fnc = @(x,t) exp(B(1)*t)*exp(1i*(x-A(1)*t));
val = zeros(length(x_nodes),1);

for i=1:100
    t = (i-1)*dt;
    for j=1:length(x_nodes)
        val(j) = fnc(x_nodes(j),t);
    end
    plot(x_nodes,val);
    xlim([0 2]);
    ylim([0.5 1]);
    xlabel('x values');
    ylabel('z(x,t)');
    title('L1');
    pause(0.5);
end


%L2
dt = 1/100;
x_nodes = linspace(0,2,100);
B = @(z)-1*z^2;
A = @(q)q;
fnc = @(x,t) exp(B(1)*t)*exp(1i*(x-A(1)*t));
val = zeros(length(x_nodes),1);

for i=1:100
    t = (i-1)*dt;
    for j=1:length(x_nodes)
        val(j) = fnc(x_nodes(j),t);
    end
    plot(x_nodes,val);
    xlim([0 2]);
    ylim([0.2 1]);
    xlabel('x values');
    ylabel('z(x,t)');
    title('L2');
    pause(0.5);
end
%}

%L3
dt = 1/100;
x_nodes = linspace(0,10,100);
B = @(z)0;
A = @(q)q+q^3;
fnc = @(x,t) exp(B(10)*t)*exp(1i*(100*x-A(100)*t));
val = zeros(length(x_nodes),1);

for i=1:100
    t = (i-1)*dt;
    for j=1:length(x_nodes)
        val(j) = fnc(x_nodes(j),t);
    end
    plot(x_nodes,val);
    xlim([0 10]);
    ylim([0 1]);
    xlabel('x values');
    ylabel('z(x,t)');
    title('L3');
    pause(0.5);
end
