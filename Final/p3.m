clear
close ALL HIDDEN
load('mimo_system')

randn('seed',1);
u1 = randn(1000,1);
randn('seed',2);
u2 = randn(1000,1);

[m, r, l] = size(MP);
u = [u1 u2];

[Y, X] =dlsim(A,B,C,D,[u1 u2]);

[Y_bar, V_bar] = YV_Form_nonzero(u,Y,10);
cap_y_hat = Y_bar*pinv(V_bar);
RSMP = recover_SYSMP(cap_y_hat,100,m,r);

plot(reshape(MP(:,1,:),[],1))
hold on
plot(reshape(RSMP(:,1,:),[],1))
title('Realized and Original Markov Parameters -- Input 1')
grid on
figure
plot(abs(reshape(MP(:,1,:)-RSMP(:,1,:),[],1)))
title('Error of Realized Markov Parameters -- Input 1')

figure
plot(reshape(MP(:,2,:),[],1))
hold on
plot(reshape(RSMP(:,2,:),[],1))
title('Realized and Original Markov Parameters -- Input 2')
grid on
figure
plot(abs(reshape(MP(:,2,:)-RSMP(:,2,:),[],1)))
title('Error of Realized Markov Parameters -- Input 2')



