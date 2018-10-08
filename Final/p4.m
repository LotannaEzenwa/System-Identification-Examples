clear
close ALL HIDDEN
load('mimo_system')

w1 = 1; % rad/sec
w2 = 2; % rad/sec
t = linspace(0,99.9,1000)';
u1 = 1 + 2*sin(w1*t) + 3*sin(w2*t);
randn('seed',2);
u2 = randn(1000,1);
u = [u1 u2];
m=2;
r=1;

[y,X] = dlsim(A,B,C,D,u);
[y_dist,X] = dlsim(A,B,C,D,u);
[Y_bar, V_bar] = YV_Form_nonzero(u2,y_dist,200);
cap_y_hat = Y_bar*pinv2(V_bar,1e-9);
RSMP = recover_SYSMP(cap_y_hat,100,m,r);

plot(reshape(MP(:,2,1:50),[],1))
hold on
plot(reshape(RSMP(:,:,1:50),[],1))
grid on
title('Realized and Original Markov Parameters -- Input 2')
figure
plot(abs(reshape(MP(:,2,1:50)-RSMP(:,:,1:50),[],1)))
title('Error of Realized Markov Parameters -- Input 2')

