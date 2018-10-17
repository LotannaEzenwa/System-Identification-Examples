clear
A_c = [0 1; -1 0];
B_c = [0;1];
C = [-1 0];
D = 1;
dt=0.2;

mat_ss = ss(A_c,B_c,C,D);
sys_dsct = c2d(mat_ss,0.2);
A_d = sys_dsct.A;
B_d = sys_dsct.B;
C_d = sys_dsct.C;
D_d = sys_dsct.D;

MarkovP = zeros(1,100);

for i=1:100
    MarkovP(i) = C_d*A_d^(i-1)*B_d;
end

clf
plot(MarkovP)
title('Lota Ezenwa -- HW1.4 -- ENGG149')
