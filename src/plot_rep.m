%DATA = csvread('Data_simulation.csv');
DATA=csvread('No_JCa_0pt18/Data_simulation.csv');
DATA2=csvread('No_JCa_0pt4/Data_simulation.csv');
close all
all_indices();
all_constants();

linestylle = ':';

a = zeros(1,43);%25 %35 Farr %43trpv
s = zeros(1,25); %27TRPV
e = zeros(1,17);
f= zeros(1,29); %24
dfdt= zeros(1,31);%24
t= zeros(3,1);
% input= zeros(1,2);

acoff   =          0;
smcoff  =          length(a);
ecoff   = smcoff + length(s);
stoff   = ecoff  + length(e);
dfdtoff = stoff   + length(f);
tijdoff = dfdtoff+ length(dfdt);
inputoff= tijdoff+ length(t);


time = DATA(:,length(DATA(1,:))-3);
time2 = DATA2(:,length(DATA2(1,:))-3);

% plots report

%min(DATA(:,acoff+flu.J_BK_k)./(DATA(:,stoff+ind.R_k)))
 total=4;
 
 
figure(1)
hold all
set(gcf,'Position', [500 62 600 700],...
        'PaperPosition', [0.634517 6.34517 20.3046 15.2284]...
        );

subplot(total,1,1)
    plot( time,DATA(:,inputoff+1)*k_C);
title('Input signal from the neuron into the synaptic cleft')
xlabel('Time [s]')
ylabel('K^+ flux [\muM m/s]')

% subplot(total,1,1)
%     plot( time,DATA(:,acoff+flu.K_s));
% title('K^+ in SC')
% xlabel('Time [s]')
% ylabel('K_s [\muM]')

subplot(total,1,2)
plot(time,DATA(:,acoff+flu.v_k)*1000)
title('Membrane Potential of the astrocyte')
xlabel('Time [s]')
ylabel('v_k [mV]')

subplot(total,1,3)
[AX2,H12,H23] = plotyy(time,DATA(:,acoff+flu.J_BK_k)./(DATA(:,stoff+ind.R_k)*VR_pa)...
    ,time,DATA(:,smcoff+flu.J_KIR_i)/(VR_ps)...
    );
xlabel('Time [s]')
set(get(AX2(1),'Ylabel'),'String','J_{BK} [\muM/s]') 
set(get(AX2(2),'Ylabel'),'String','J_{KIR} [\muM/s]')
%legend('BK-channel','KIR-channel')
title('The contribution of the BK- and KIR-channel to K_p')

subplot(total,1,4)
plot(time,DATA(:,stoff+ind.K_p)...
     )
title('Potassium concentration in perivascular space')
xlabel('Time [s]')
ylabel('K_p  [\muM]')

%% DFDT plot

figure(4)
set(gcf,'Name','DFDT')
set(gcf,'Position', [24 62 1616 904],...
        'PaperPosition', [0.634517 6.34517 20.3046 15.2284]...
        );
subplot(3,3,1)
hold on
plot(time, DATA(:,stoff+ind.Ca_i) )
plot(time2, DATA2(:,stoff+ind.Ca_i),'r' )
h6(1) = gca();
xlabel('time in s')
ylabel('Ca_i in uM')
title('SMC [Ca^{2+}]')
hold all

subplot(3,3,2)
hold on
plot(time, DATA(:,stoff+ind.Ca_j) )
plot(time2, DATA2(:,stoff+ind.Ca_j),'r' )
h6(2) = gca();
xlabel('time in s')
ylabel('Ca_j in uM')
title('EC [Ca^{2+}]')
hold all

subplot(3,3,3)
hold on
plot(time, DATA(:,stoff+ind.I_i) )
plot(time2, DATA2(:,stoff+ind.I_i),'r' )
h6(8) = gca();
xlabel('time in s')
ylabel('I_i in uM')
title('SMC [IP_{3}]')
hold all

subplot(3,3,4)
hold on
plot(time, DATA(:,stoff+ind.I_j) )
plot(time2, DATA2(:,stoff+ind.I_j),'r' )
h6(9) = gca();
xlabel('time in s')
ylabel('I_j in uM')
title('EC [IP_{3}]')
hold all

subplot(3,3,5)
hold on
plot(time, DATA(:,stoff+ind.s_i) )
plot(time2, DATA2(:,stoff+ind.s_i),'r' )
h6(3) = gca();
xlabel('time in s')
ylabel('s_i in uM')
title('SR [Ca^{2+}]')
hold all

subplot(3,3,6)
hold on
plot(time, DATA(:,stoff+ind.s_j) )
plot(time2, DATA2(:,stoff+ind.s_j),'r' )
h6(4) = gca();
xlabel('time in s')
ylabel('s_j in uM')
title('ER [Ca^{2+}]')
hold all

subplot(3,3,7)
hold on
plot(time, DATA(:,stoff+ind.v_i) )
plot(time2, DATA2(:,stoff+ind.v_i),'r' )
h6(5) = gca();
xlabel('time in s')
ylabel('v_i in mV')
title('SMC membrane voltage')
hold all

subplot(3,3,8)
hold on
plot(time, DATA(:,stoff+ind.v_j) )
plot(time2, DATA2(:,stoff+ind.v_j),'r' )
h6(6) = gca();
xlabel('time in s')
ylabel('v_j in mV')
title('EC membrane voltage')
hold all

subplot(3,3,9)
hold on
plot(time, DATA(:,stoff+ind.w_i) )
plot(time2, DATA2(:,stoff+ind.w_i),'r' )
h6(7) = gca();
xlabel('time in s')
ylabel('w_i [-]')
title('open probability K^+ channel')
hold all


linkaxes(h6, 'x');


figure(9)
hold all
set(gcf,'Name','TRPV results')
set(gcf,'Position', [500 62 600 700],...
        'PaperPosition', [0.634517 6.34517 20.3046 15.2284]...
        );
subplot(2,1,1)
hold on
plot(time, DATA(:,stoff+ind.Ca_p) )
plot(time2, DATA2(:,stoff+ind.Ca_p),'r' )
xlabel('time in s')
ylabel('Ca_i in uM')
title('PVS [Ca^{2+}]')
subplot(2,1,2)
hold on
    plot( time,DATA(:,flu.J_TRPV_k));
    plot( time2,DATA2(:,flu.J_TRPV_k),'r');
title('Flux through TRPV-channel')
xlabel('Time [s]')
ylabel('J_TRPV_k [\muM ps]')

% subplot(3,1,3)
%     plot( time,DATA(:,smcoff+flu.J_Ca));
% title('Flux through Ca-channel')
% xlabel('Time [s]')
% ylabel('J_Ca [\muM ps]')



%% Potassium plots

% figure;
% subplot(2,2,1)
% plot(time,state(:,ind.K_i))
% title('K^+ in SMC')
% subplot(2,2,2)
% plot(time,DATA(:,smcoff+flu.J_KIR_i))
% title('KIR channel outflux')
% subplot(2,2,3)
% plot(time,-DATA(:,smcoff+flu.J_K_i))
% title('J K_i outflux')
% subplot(2,2,4)
% plot(time,DATA(:,smcoff+flu.J_NaK_i))
% title('J NaK_i influx')