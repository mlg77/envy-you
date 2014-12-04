clean
tic

% global variables
global CASE J_PLC startpulse lengthpulse C_Hillmann stretch_ch only_Koenig NVU Glu_start Glu_end wss_start wss_end

%% NO pathway
global m %(cGMP coupling (0 - lowest influence to 2 - highest influence))
m = 2;
% lalaa = 1;
% 
% cai = [];
% logcai = [];
% Po = [];
% vca3 = [];
% rcgmp1 = [];
% vii = [];
% timee = [];

% for vivi = -100:10:200
%% Parameters to adjust the model:
t_start = 0;
t_end = 700;
startpulse  = 20000;  % (s) 
lengthpulse = 200000;  % (s) 
Glu_start   = 200;
Glu_end     = 400;
wss_start   = 100000; 
wss_end     = 120000;
CASE        = 2;    % (see all_constants.m for details)
J_PLC 		= 0.18;  % 0.18(steady) %0.4(fluctuating) (muM s-1) EC agonist concentration  
C_Hillmann  = 1;    % scaling factor for the Hai&Murphy rate constants (see all_constants.m for details)
stretch_ch  = 'ON'; % choose 'ON'/'OFF' to activate/deactivate stretch-activated channels in EC and SMC
only_Koenig = 'OFF';% choose 'ON'/'OFF' to simulate only the Koenigsberger model (other sub-models will still be considered, but the KIR channel is set to 0)
NVU         = 2;     % 1=NVU 1.0 , 2=NVU 1.1, 3=NVU 1.0 + EET, 4= NVU 1.0 + Ca2+
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
options = odeset('OutputFcn',@odeprogWD,'Events',@odeabort,'Stats','on','RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1); 
[t,state] = ode15s(@DEsyst,[t_start t_end],state0,options);

%% Write output and info to file/cmd
output.info.completiontime = toc;
fprintf('ODE solution time: %.3f seconds\n', output.info.completiontime)

%% Plot statement:
plot_all()
hold all

%% save figures & parameters
%save_all()


% to create .tikz figures:
% matlab2tikz('test.tikz', 'height', '\figureheight', 'width', '\figurewidth');

% 
% % figure
% % subplot(4,1,1)
% % plot(time,DATA(:,smcoff+flu.Kactivation_i),'-x')
% % subplot(4,1,2)
% % plot(time,DATA(:,smcoff+flu.v_Ca3),'-x')
% % subplot(4,1,3)
% % plot(time,state(:,ind.Ca_i),'-x')
% % subplot(4,1,4)
% % plot(time,DATA(:,smcoff+flu.R_cGMP1),'-x')
% 
% foo0 = state(:,ind.Ca_i);
% foo1 = DATA(:,smcoff+flu.Kactivation_i);
% foo2 = DATA(:,smcoff+flu.v_Ca3);
% foo3 = DATA(:,smcoff+flu.R_cGMP1);
% foo4 = state(:,ind.v_i);
% 
% cai{lalaa} = mean(foo0(200:end));
% % logcai{lalaa} = log10(cai);
% Po{lalaa} = mean(foo1(200:end));
% vca3{lalaa} = mean(foo2(200:end));
% rcgmp1{lalaa} = mean(foo3(200:end));
% vii{lalaa} = mean(foo4(200:end));
% timee{lalaa} = time;
% 
% clear foo0 foo1 foo2 foo3 foo4
% 
% lalaa = lalaa + 1
% end
% % vi(lalaa)   = vivi
% % Po(lalaa)   = DATA(end,smcoff+flu.Kactivation_i)
% % vca3(lalaa) = DATA(end,smcoff+flu.v_Ca3)
% % cai(lalaa)  = state(end,ind.Ca_i) 
% % 
% % 
% % lalaa = lalaa+1
% % end

% figure; plot(time,DATA(:,smcoff+flu.M)+state(:,ind.AMp)+state(:,ind.AM)+state(:,ind.Mp)); hold on;
% plot(time,DATA(:,smcoff+flu.M),'r'); plot(time,state(:,ind.Mp),'g'); plot(time,state(:,ind.AMp),'b');plot(time,state(:,ind.AM),'k');
% legend('Total Myosin','[M]','[Mp]','[AMp]','[AM]')
% title('New')

% to plot a single flux, type in plot(time,DATA(:,flu.(name))     
% to plot a single state variable, type in plot(time,DATA(:,ind.(name))
%(don't forget to put the offset!! e.g. smcoff+flu.1_c)
