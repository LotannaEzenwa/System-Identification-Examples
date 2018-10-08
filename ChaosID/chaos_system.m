clear;

t = 0:0.01:10;
q = 2;
a = sin(q*t);

B = [
    0 0
    1 0 
    0 1
    ];
C = [
    0 1 0
    0 0 1
    ];
D = 0;

lp = @(d) ([d; sin(d)]);
plot(t,lp(t));

