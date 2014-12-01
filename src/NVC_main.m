clean
tic
global Ca_0 v_0 K_0 R_0 Rk_0 T_c nonDim
nonDim      = 'ON'; % choose 'ON'/'OFF' to activate/deactivate non-dimensionalisation of the model

if strcmp(nonDim,'ON') == 1
    %% This is the non-dimensionalised version of the envy-you model - version 1.1
    % Characteristic dimensional quantities are:
    Ca_0 = 0.8; % microM  (Ca2+ concentration)
    v_0  = -40; % mV      (membrane potential)
    K_0  = 1.4e4;% microM (K+ concentration)
    R_0  = 2.6e-5; % m    (Radius)
    Rk_0 = 0.7e-7; % m    (AC volume-area ratio)
elseif strcmp(nonDim,'OFF') == 1  
    Ca_0 = 1;
    v_0  = 1;
    K_0  = 1;
    R_0  = 1;
    Rk_0 = 1;
end

%% load the constants for the fluxes and pointers:
all_indices();
all_constants();

%% characteristic time scale:
T_c  = Ca_0/C_i; % = 0.0145 s  - can be changed!

%% global variables
global CASE J_PLC startpulse lengthpulse C_Hillmann stretch_ch only_Koenig NVU

%% Parameters to adjust the model:
t_start = 0/T_c;
t_end = 400/T_c;
startpulse  = 200/T_c;  % (s) 
lengthpulse = 100/T_c;  % (s) 
CASE        = 2;    % (see all_constants.m for details)
J_PLC 		= 0.18;  % 0.18(steady) %0.4(fluctuating) (muM s-1) EC agonist concentration  
C_Hillmann  = 1;    % scaling factor for the Hai&Murphy rate constants (see all_constants.m for details)
stretch_ch  = 'ON'; % choose 'ON'/'OFF' to activate/deactivate stretch-activated channels in EC and SMC
only_Koenig = 'OFF';% choose 'ON'/'OFF' to simulate only the Koenigsberger model (other sub-models will still be considered, but the KIR channel is set to 0)
NVU         = 1;     % 1=NVU 1.0 , 2=NVU 1.1, 3=NVU 1.0 + EET, 4= NVU 1.0 + Ca2+

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
% plot_all()
% hold all

%% save figures & parameters
%save_all()
%%state:
for i = 1:24
   figure(ceil(i/6));
   spl_no = mod(i,6);
   if spl_no == 0
      spl_no = 6;
   end
   subplot(2,3,spl_no)
   plot(t,state(:,i))
   ylabel(i)
end

%% to create .tikz figures:
% matlab2tikz('test.tikz', 'height', '\figureheight', 'width', '\figurewidth');


% figure; plot(time,DATA(:,smcoff+flu.M)+state(:,ind.AMp)+state(:,ind.AM)+state(:,ind.Mp)); hold on;
% plot(time,DATA(:,smcoff+flu.M),'r'); plot(time,state(:,ind.Mp),'g'); plot(time,state(:,ind.AMp),'b');plot(time,state(:,ind.AM),'k');
% legend('Total Myosin','[M]','[Mp]','[AMp]','[AM]')
% title('New')

% to plot a single flux: plot(time,DATA(:,flu.(name))  - (don't forget the offset!! e.g. smcoff+flu.1_c)      
% to plot a single state variable: plot(time,state(:,ind.(name))



%%  plot getRef input:
% startpulse = 20;
% lengthpulse = 50;
% 
% for i = 1:100
%     t(i) = i;
%     y(i) = getRef(i,'ft');
% end
% figure; plot(t,y)
% 
% for i = 1:100
%     t(i) = i;
%     y(i) = getRef(i,'fluxft');
% end
% figure; plot(t,y)