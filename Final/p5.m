p2;
close ALL HIDDEN
load('mimo_system')
ts = 0.1;

[A_s,B_s,C_s,D_s] = d2cm(A_r,B_r,C_r,D,ts,'zoh');


Q_star = rand(size(C_s));
Q = [C_s;Q_star];
[w,q] = size(A_s);

A_star = Q*A_s*inv(Q);
B_star = Q*B_s;
C_star = C_s*inv(Q);

n = w/2;

T11 = eye(n);
T12 = zeros(n);
T21 = A_star(1:n,1:n);
T22 = A_star(1:n,n+1:end);



T = [
    T11 T12
    T21 T22
    ];

A_pr = T*A_star*inv(T);
B_pr = T*B_star;
C_pr = C_star*inv(T);

[n, n2] = size(C_pr);


X = -A_pr(n+1:end,1:n);
Y = -A_pr(n+1:end,n+1:end);
Z = B_pr(n+1:end,:);
I_n = eye(n);
I_nn = eye(n*n);
O_n = zeros(n);
O_nn = zeros(n*n);
O_2n1 = zeros(n*n,1);

R = [
    kron(X',I_n) O_nn -I_nn
    kron(Y',I_n) -I_nn O_nn 
    kron(Z',I_n) O_nn O_nn
    kron(X',I_n)-kron(I_n,X') O_nn O_nn
    kron(Y',I_n)-kron(I_n,Y') O_nn O_nn
    ];



G = [
    O_2n1
    O_2n1
    Z(:)
    O_2n1
    O_2n1
    ];

P = pinv(R)*G;
eqm = reshape(P,2,2,3);

