clean
tic

% global variables
global CASE J_PLC startpulse lengthpulse C_Hillmann stretch_ch only_Koenig

%% Parameters to adjust the model:
t_start     = 0;
t_end       = 200;
startpulse  = 1000;  % (s) 
lengthpulse = 40;  % (s) 
CASE        = 2;    % (see all_constants.m for details)
J_PLC 		= 0.18;  % (muM s-1) EC agonist concentration  
C_Hillmann  = 1;    % scaling factor for the Hai&Murphy rate constants (see all_constants.m for details)
stretch_ch  = 'ON'; % choose 'ON'/'OFF' to activate/deactivate stretch-activated channels in EC and SMC
only_Koenig = 'OFF';% choose 'ON'/'OFF' to simulate only the Koenigsberger model (other sub-models will still be considered, but the KIR channel is set to 0)

%% load the constants for the fluxes and pointers:
all_indices();
all_constants();
%% load the initial conditions of the system:
state0 = InitCond();
%% Ensure single filenames for the writing of data in other files
global csvfilename
csvfilename = 'Data_simulation.csv';
try
delete(csvfilename) % remove file, if present from older simulation.
end
%% Solve the proces from initial position tot Steady State:
options = odeset('OutputFcn',@odeprogWD,'Events',@odeabort,'Stats','on','RelTol', 1e-3, 'AbsTol', 1e-03, 'MaxStep', 1); 
[t,state] = ode15s(@DEsyst,[t_start t_end],state0,options);

%% Write output and info to file/cmd
output.info.completiontime = toc;
fprintf('ODE solution time: %.3f seconds\n', output.info.completiontime)

%% Plot statement:
plot_all()
hold all

%% plot all of neuron model's state variables:
figure(8);
set(gcf,'Name','Neuron Membrane potentials')
for i = 25:26
subplot(1,2,i-24)
plot(time,state(:,i))
ylabel(i)
end

figure(9);
set(gcf,'Name','Neuron Gating variables')
for i = 27:38
subplot(4,4,i-26)
plot(time,state(:,i))
ylabel(i)
end

figure(10);
set(gcf,'Name','Neuron Concentration of ions')
for i = 39:49
subplot(3,4,i-38)
plot(time,state(:,i))
ylabel(i)
end

%% plot all of neurons model's fluxes:

figure(11); 
set(gcf,'Name','Neuron Fluxes1-20')
for i = 1:20
subplot(4,5,i); 
plot(time,DATA(:,neoff+i)); 
ylabel(i)
end

figure(12); 
set(gcf,'Name','Neuron Fluxes21-40')
for i = 21:40
subplot(4,5,i-20); 
plot(time,DATA(:,neoff+i)); 
ylabel(i)
end

figure(13); 
set(gcf,'Name','Neuron Fluxes41-63')
for i = 41:62
subplot(5,5,i-40); 
plot(time,DATA(:,neoff+i)); 
ylabel(i)
end
figure(14); 
set(gcf,'Name','Neuron Fluxes64-66')
for i =63:66
subplot(3,2,i-62); 
plot(time,DATA(:,neoff+i)); 
ylabel(i)
end
%% Extra plots of neuron model
% figure(15)
% set(gcf,'Name','Membrane potential ')
% plot(time, state(:,ind.N_sa_Em) )
% xlabel('time in s')
% ylabel('Soma potential in mV')
% title('Membrane potential')
% 
% figure(16)
% set(gcf,'Name','Extrcellular potassium and sodium current ')
% subplot(2,1,1);
% plot(time, state(:,ind.N_e_K) )
% xlabel('time in s')
% ylabel('ECS K+ in mM')
% title('Extracellular potassium')
% subplot(2,1,2);
% plot(time, state(:,ind.N_e_Na) )
% xlabel('time in s')
% ylabel('ECS Na+ in mM')
% title('Extracellular Sodium')
% 
% figure(17)
% set(gcf,'Name','Neuron_to_Radius')
% subplot(3,1,1)
% plot( time,0.001*DATA(:,acoff+flu.K_s));
% title('[K^+] in synaptic cleft')
% xlabel('Time [s]')
% ylabel('[K^+]_s [mM]')
% hold all
% subplot(3,1,2)
% plot(time, state(:,ind.Ca_i) )
% title('[Ca^{2+}] in smooth muscle cell')
% xlabel('Time [s]')
% ylabel('[Ca^{2+}]_i [\muM]')
% hold all
% 
% 
% 
% subplot(3,1,3)
% plot(time,1e6*state(:,ind.R))
% title('Radius')
% xlabel('Time [s]')
% ylabel('Radius [\mum]')
% %axis([0 250 15 27])
% hold all
% % axis([0 250 15 27])
% 
% 
% 
% figure(18)
% set(gcf,'Name','Soma Membrane potential')
% plot(time,state(:,25))
% xlabel('time in s')
% ylabel('Membrane potential in mV')
% title('Soma Membrane potential ')
% 
% figure(19)
% set(gcf,'Name','Extrcellular potassium ')
% plot(time, state(:,ind.N_e_K) )
% xlabel('time in s')
% ylabel('ECS K^+ in mM')
% title('Extracellular potassium')
% 
% figure(20)
% set(gcf,'Name','[Ca^{2+}] in smooth muscle cell ')
% plot(time, state(:,ind.Ca_i) )
% title('[Ca^{2+}] in smooth muscle cell')
% xlabel('Time [s]')
% ylabel('[Ca^{2+}]_i [\muM]')
% 
% figure(21)
% set(gcf,'Name','Radius of the artery ')
% plot(time,1e6*state(:,ind.R))
% title('Radius of the artery')
% xlabel('Time [s]')
% ylabel('Radius [\mum]')
% 
% figure(22)
% set(gcf,'Name','Cerebral Blood Flow ')
% plot(time, DATA(:,neoff+flu.V_CBF) )
% xlabel('time in s')
% ylabel('CBF in mM/s')
% title('Cerebral Blood Flow')
% 
% figure(23)
% set(gcf,'Name','Tissue oxygen Concentration ')
% plot(time, state(:,ind.N_O2) )
% xlabel('time in s')
% ylabel('Tissue O2 in mM')
% title('Tissue oxygen Concentration ')

%% save figures & parameters
save_all()


% to create .tikz figures:
% matlab2tikz('test.tikz', 'height', '\figureheight', 'width', '\figurewidth');


% figure; plot(time,DATA(:,smcoff+flu.M)+state(:,ind.AMp)+state(:,ind.AM)+state(:,ind.Mp)); hold on;
% plot(time,DATA(:,smcoff+flu.M),'r'); plot(time,state(:,ind.Mp),'g'); plot(time,state(:,ind.AMp),'b');plot(time,state(:,ind.AM),'k');
% legend('Total Myosin','[M]','[Mp]','[AMp]','[AM]')
% title('New')

% to plot a single flux, type in plot(time,DATA(:,flu.(name))     
% to plot a single state variable, type in plot(time,state(:,ind.(name))
%(don't forget to put the offset!! e.g. smcoff+flu.1_c)
