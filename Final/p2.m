clear
close ALL HIDDEN
load('mimo_system')


n = 45;
[m,r,l] = size(MP);

desired=4;

[H_0,H_1] = mimo_hankel(MP,desired/m);

[U,S,V] = svd(H_0);
sing_n = diag(S((S./max(max(S)))>1e-13));

n_e = length(nonzeros(sing_n));

Un = U(:, 1:n_e);
Vn = V(:, 1:n_e);


A_r = sing_n^-0.5*Un'*H_1*Vn*sing_n^-0.5;
B_0 = sing_n^0.5*Vn';
B_r = B_0(:,1:r);
C_0 = Un*sing_n^0.5;
C_r = C_0(1:m,:);

era_MP = zeros(size(MP));

for i=2:l
    era_MP(:,:,i) = C_r*A_r^(i-2)*B_r;
end

hankel_matrix_size = size(H_0)
semilogy(diag(sing_n./(max(max(sing_n)))),'-*','color','c')
title('Normalized Singular Values')
figure

plot(era_MP(:));
hold on
plot(MP(:),'color','g');
title('Realized and Original Markov Parameters')
xlabel('Stacked Markov Parameter #')
figure
plot(abs(era_MP(:)-MP(:)))
title('Error of Realized Markov Parameters')


