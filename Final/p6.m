clear
load('hubbledata.mat')
p=200;

[Y_bar, V_bar] = YV_Form_nonzero(u,10*y,p);
cap_y_hat = Y_bar*pinv2(V_bar,1e-15);

RSMP = recover_SYSMP(cap_y_hat,100,q,r);
desired = q*p;


[H_0,H_1] = mimo_hankel(RSMP,desired/q);
H_0 = H_0(1:desired,1:desired);
H_1 = H_1(1:desired,1:desired);

[U,S,V] = svd(H_0);
sing_n = diag(S((S./max(max(S)))>1e-13));

n_e = length(nonzeros(sing_n));

Un = U(:, 1:n_e);
Vn = V(:, 1:n_e);


A_r = sing_n^-0.5*Un'*H_1*Vn*sing_n^-0.5;
B_0 = sing_n^0.5*Vn';
B_r = B_0(:,1:r);
C_0 = Un*sing_n^0.5;
C_r = C_0(1:q,:);

era_MP = zeros(size(RSMP));
l = 50

for i=2:l
    era_MP(:,:,i) = C_r*A_r^(i-2)*B_r;
end

hankel_matrix_size = size(H_0)
semilogy(diag(sing_n./(max(max(sing_n)))),'-*')
title('Normalized Singular Values')

[Y_2,X_2] = dlsim(A_r,B_r,C_r,0,u);
subplot(3,1,1);
plot(y(:,1)-Y_2(:,1))
subplot(3,1,2);
plot(y(:,2)-Y_2(:,2))
subplot(3,1,3);
plot(y(:,3)-Y_2(:,3))



