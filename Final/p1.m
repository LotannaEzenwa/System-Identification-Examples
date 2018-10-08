%% ENGS149 Final Exam
%% Problem 1

clear
m1 = 1;
m2 = 1;
c1 = 3;
c2 = 2;
k1 = 30;
k2 = 20;

A_c = [
    0 0 1 0
    0 0 0 1
    -(k1+k2)/m1 k2/m1 -(c1+c2)/m1 c2/m1
    k2/m2 -k2/m2 c2/m2 -c2/m2
    ];
B_c = [
    0 0
    0 0
    1/m1 0  
    0 1/m2
    ];

C_c = [
    1 0 0 0
    0 1 0 0
    ];
D_c = [
    0 0 
    0 0
    ];

[n_st, l] = size(A_c);
[o,n_in] = size(B_c);
[n_out,p] = size(C_c);

dt = 0.1; % seconds

c_system = ss(A_c,B_c,C_c,D_c);

d_system = c2d(c_system,dt);
A = d_system.A;
B = d_system.B;
C = d_system.C;
D = d_system.D;

n = 100;

x = zeros(n_st,n);
x_b1 = x;
x_b2 = x;
u = zeros(n_in,n);
u_b1 = u;
u_b2 = u;
u_b1(:,1) = [1;0];
u_b2(:,1) = [0;1];
y_1 = zeros(2,n);
y_2 = zeros(2,n);
MP = zeros(n_in,n_out,n);

MARKOV_P = zeros(n_in,n_out,n);


MARKOV_P(:,:,1) = D;

for i = 1:n
    x_b1(:,i+1) = A*x_b1(:,i) + B*u_b1(:,i);
    x_b2(:,i+1) = A*x_b2(:,i) + B*u_b2(:,i);
    y_1(:,i) = C*x_b1(:,i) + D*u_b1(:,i);
    y_2(:,i) = C*x_b2(:,i) + D*u_b2(:,i);
    MP(:,1,i) = y_1(:,i);
    MP(:,2,i) = y_2(:,i);
    MARKOV_P(:,:,i+1) = C*A^(i-1)*B;
end

subplot(2,2,1)
plot(reshape(MARKOV_P(1,1,:),[],1))
hold on
grid on
plot(reshape(MP(1,1,:),[],1),'*')
title('Markov Parameters vs. Impulse Response, Input 1, Output 1')
subplot(2,2,2)
plot(reshape(MARKOV_P(1,2,:),[],1))
hold on
grid on
plot(reshape(MP(1,2,:),[],1),'*')
title('Markov Parameters vs. Impulse Response, Input 2, Output 1')

subplot(2,2,3)
plot(reshape(MARKOV_P(2,1,:),[],1))
hold on
grid on
plot(reshape(MP(2,1,:),[],1),'*')
title('Markov Parameters vs. Impulse Response, Input 1, Output 2')

subplot(2,2,4)
plot(reshape(MARKOV_P(2,2,:),[],1))
hold on
grid on
plot(reshape(MP(2,2,:),[],1),'*')
title('Markov Parameters vs. Impulse Response, Input 2, Output 2')


save('mimo_system',"A","B","C","D","MP","d_system",'A_c','B_c','C_c','D_c')

