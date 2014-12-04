K_mO2_j = 7.7; % uM Chen2006
K_mO2_n = 243; % 140-145 uM Chen2007

for i = 1:200;
    oxygen(i) = i;
    P_NO_j(i) = 0.055 * oxygen(i) / (K_mO2_j + oxygen(i));
    P_NO_n(i) = 3.8 * oxygen(i) / (K_mO2_n + oxygen(i));
end

figure; plot(oxygen,P_NO_j)
xlabel('[O_2] in EC [uM]')
ylabel('NO production rate in EC [uM/s]')
figure; plot(oxygen,P_NO_n)
xlabel('[O_2] in NC [uM]')
ylabel('NO production rate in NE [uM/s]')
