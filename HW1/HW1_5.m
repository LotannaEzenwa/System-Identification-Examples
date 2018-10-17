HW1_3
close ALL HIDDEN
T = [1 2; 1 0];

MarkovPT = zeros(1,100);
TA = T*A_d*inv(T);
TB = T*B_d;
TC = C_d*inv(T);

for i=1:100
    MarkovPT(i) = TC*TA^(i-1)*TB;
end

clf
hold on
plot(abs(MarkovPT - MarkovP))
title('Lota Ezenwa -- HW1.5 -- ENGG149')
legend('MP_t - MP')
