close all;

syms x;
y1 = triangularPulse(-1.5,-0.5,x);
y2 = triangularPulse(-1,0,x);
y3 = triangularPulse(-0.5,0.5,x);
y4 = triangularPulse(0,1,x);
y5 = triangularPulse(0.5,1.5,x);

tn=-1:0.01:1;
y1n = double( subs(y1, x, tn)) ;
y2n = double( subs(y2, x, tn)) ;
y3n = double( subs(y3, x, tn)) ;
y4n = double( subs(y4, x, tn)) ;
y5n = double( subs(y5, x, tn)) ;

figure
plot(tn,y1n,tn,y2n,tn,y3n,tn,y4n,tn,y5n);

syms x;
y1 = chebyshevT(0,x);
y2 = chebyshevT(1,x);
y3 = chebyshevT(2,x);
y4 = chebyshevT(3,x);
y5 = chebyshevT(4,x);

tn=-1:0.01:1;
y1n = double( subs(y1, x, tn)) ;
y2n = double( subs(y2, x, tn)) ;
y3n = double( subs(y3, x, tn)) ;
y4n = double( subs(y4, x, tn)) ;
y5n = double( subs(y5, x, tn)) ;

figure
plot(tn,y1n,tn,y2n,tn,y3n,tn,y4n,tn,y5n);
