%% General constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Farad       = 96500         ;% [C mol-1]      Faradays constant
R_gas       = 8.315         ;% [J mol-1K-1]
Temp        = 300           ;% [K]

unitcon    = 10^(3)         ;% [-]            Factor to convert equations to another unit

%% Constants of Ostby, the values come from CELL ML
 
L_p         = 2.1e-9        ;% [m uM-1s-1]
R_tot       = 8.79e-8       ;% [m]
X_k         = 12.41e-3      ;% [uMm]
z_Na        = 1             ;% [-]
z_K         = 1             ;% [-]
z_Cl        = -1            ;% [-]
z_NBC       = -1            ;% [-]
g_K_k       = 40            ;% [ohm-1m-2]
g_KCC1_k    = 1e-2          ;% [ohm-1m-2]
g_NBC_k     = 7.57e-1       ;% [ohm-1m-2]
g_Cl_k      = 8.797e-1      ;% [ohm-1m-2]
g_NKCC1_k   = 5.54e-2       ;% [ohm-1m-2]
g_Na_k      = 1.314         ;% [ohm-1m-2]
J_NaK_max   = 1.42e-3       ;% [uMm s-1]
K_Na_k      = 10e3          ;% [uM]
K_K_s       = 1.5e3         ;% [uM]
k_C         = 7.35e-5       ;% [muM s-1]

% Constants of Neuron Model

Ra          = 1.83e5     ;% [ohms] Input reisitance of dendritic tree 
dhod        = 4.5e-2     ;% [cm] Half length of dendrite
As          = 1.586e-5   ;% [cm2] Surface area of soma 
Ad          = 2.6732e-4  ;% [cm2] Surface area of dendrite
Vs          = 2.160e-9   ;% [cm3] Volume of soma
Vd          = 5.614e-9   ;% [cm3] Volume of dendrite
Cm          = 7.5e-5     ;% [s/Ohmcm2] Membrane capacitance
Imax        = 1.48e-3    ;% [mA/cm2] Na+/K+ -ATPase rate
R           = 8.31       ;% [mV coulomb/mmol K] Universal gas constant
T           = 310        ;% [K] Absolute temperature
F           = 96.485     ;% [coulomn/mmol] Faraday constant
ph          = 26.6995    ;% R*T/F
fe1         = (0.15*Vs)/(Vs+Vd)  ;% ve/(Vs+Vd)
fe2         = (0.15*Vd)/(Vs+Vd) ;
gNaP_GHk    = 0.5e-6       ;% [mA cm]
gKDR_GHk    = 2.5e-5      ;% [mA cm]
gKA_GHk     = 0.25e-5       ;% [mA cm]
gNMDA_GHk   = 0.25e-5       ;% [mA cm]
% gNaP_GHk    = 2e-6       ;% [mA cm]
% gKDR_GHk    = 10e-5      ;% [mA cm]
% gKA_GHk     = 1e-5       ;% [mA cm]
% gNMDA_GHk   = 1e-5       ;% [mA cm]
gNaleak_sa  = 9.6605e-7  ;
gKleak_sa   = 3.1292e-6  ;
gClleak_sa  = 10*gNaleak_sa;
gNaleak_d   = 9.6064e-7  ;
gKleak_d    = 3.1279e-6 ;
gClleak_d   = 10*gNaleak_d;
O2_0          = 2e-2       ;% [mM]
alph        = 0.05       ;% Percentage of ATP production that is independent of Oxygen
D_Na        = 1.33e-5   ;% cm^2/s
D_K         = 1.96e-5   ;% cm^2/s
D_Cl        = 2.03e-5   ;% cm^2/s
N_e_K_init  = 3.5        ;% mM concentration of K in the extracellular space
N_sa_Na_init= 10         ;% mM concentration of Na inside the soma
N_d_Na_init = 10         ;% mM concentration of Na inside the dendrite
Mu          = 8e-6          ;%ms-1
B0          = 200           ;%mM Effective total buffer concentration

%%Oxygen Depndency

LU_R_init  =15e-6;
CBF_init      =2.5e-2;
O2_b         =4e-2;
gam =0.1;

%% Constants of the BK-channel

A_ef_k      = 3.7e-9        ;% m2       Area of an endfoot of an astrocyte, equal to Area astrocyte at synaptic cleft
v_6         = 22e-3         ;% V        Constant
v_4         = 14.5e-3       ;% V        A measure of the spread of the distrubution
psi_w       = 2.664         ;% s-1      A characteristic time
G_BK_k      = 4.3e3         ; % pS      Constant estimation based on Ermentrout
g_BK_k      = G_BK_k*10^(-12)/A_ef_k ;% ohm-1m-2  Specific capacitance of the BK-Channel in units of Ostby
VR_pa       = 0.001       ;% [-]       The estimated volume ratio of perivascular space to astrocyte: Model estimation
VR_ps       = 0.001       ;% [-]       The estimated volume ratio of perivascular space to SMC: Model Estimation	




%% SMC constants


F_il = 7.5e2            ;%[-] scalingsfactor to fit the experimental data of Filosa
z_1 =4.5                ;%[-] parameter fitted on experimental data of Filosa
z_2 =-1.12e2            ;%[-] parameter fitted on experimental data of Filosa
z_3 =4.2e-1             ;%[-] parameter fitted on experimental data of Filosa
z_4 =-1.26e1            ;%[-] parameter fitted on experimental data of Filosa
z_5 =-7.4e-2            ;%[-] parameter fitted on experimental data of Filosa

% mVmicroM-1 The change in membrane potential by a scaling factor

% Koeningsberger et al.

Fmax_i		= 0.23;		% [microM/s]
Kr_i 		= 1; 		% [microM] ; Half saturation constant for agonist-dependent Ca$^{2+}$ entry
G_Ca		= 0.00129;	% [microM/mV/s]
v_Ca1		= 100;		% [mV]
v_Ca2		= -24;	    % [mV]  why did we change it to -35 temporarely???
R_Ca		= 8.5;		% [mV]
G_NaCa		= 0.00316;	% microM/mV/s
c_NaCa		= 0.5;		% microM
v_NaCa		= -30;
B_i			= 2.025;
cb_i		= 1;
C_i			= 55;
sc_i		= 2;
cc_i		= 0.9;
D_i			= 0.24;
vd_i		= -100;
Rd_i		= 250;
L_i			= 0.025;
gam			= 1970; % mVmicroM-1 The change in membrane potential by a scaling factor
F_NaK		= 0.0432;
G_Cl		= 0.00134;
v_Cl		= -25;
G_K			= 0.00446;
vK_i		= -94; 
lab 		= 45;
c_w			= 0;
bet			= 0.13;
v_Ca3		= -27;
R_K			= 12;
k_i			= 0.1;
K_d         = 1;            % = 1e3 nM Gonzalez
B_T         = 100;          % = 1e5 nM Gonzalez
Kinf_i      = 1e5;          % 100 mM K+ concentration in SMC

G_stretch   = 0.0061;       % uM mV-1 s-1
P_str       = 30;
Esac        = -18;          % mV
alpha1      = 0.0074;
sig0        = 500;


%% EC constants
% EC Koeningsberger et al.

Fmax_j		= 0.23;		% [microM/s]
Kr_j		= 1;
B_j 		= 0.5;
cb_j		= 1;
C_j			= 5;
sc_j		= 2;
cc_j		= 0.9;
D_j			= 0.24;
L_j			= 0.025;
G_cat 		= 0.66e-3;
E_Ca		= 50;
m3cat		= -0.18; %-6.18; %changed value!!! 
m4cat 		= 0.37;
J0_j 		= 0.029; %constant Ca influx (EC)
C_m 		= 25.8;
G_tot		= 6927;
vK_j 		= -80;
a1			= 53.3;
a2			= 53.3;
b			= -80.8;
c 			= -0.4; %-6.4; %changed value!!! 
m3b			= 1.32e-3;
m4b			= 0.3;
m3s			= -0.28;
m4s			= 0.389;
G_R			= 955;
v_rest		= -31.1;
k_j			= 0.1;





global CASE J_PLC g_hat p_hat p_hatIP3

if CASE==0
g_hat 		= 0;
p_hat 		= 0;
p_hatIP3 	= 0;
elseif CASE==1
g_hat 		= 0.5;
p_hat 		= 0;
p_hatIP3 	= 0.05;
elseif CASE==2
g_hat 		= 0.5;
p_hat 		= 0.05;
p_hatIP3 	= 0.05;
elseif CASE==3
g_hat 		= 0;
p_hat 		= 0;
p_hatIP3 	= 0.05;
elseif CASE==4 % was CASE 5 before!
g_hat 		= 0.5;
p_hat 		= 0.05;
p_hatIP3 	= 0;
elseif CASE==5 % was CASE 6 before!
g_hat 		= 0.5;
p_hat 		= 0;
p_hatIP3 	= 0;
elseif CASE==6 % was CASE 7 before!
g_hat 		= 0;
p_hat 		= 0.05;
p_hatIP3 	= 0;
elseif CASE==7 % was CASE 8 before!
g_hat 		= 0;
p_hat 		= 0.05;
p_hatIP3 	= 0.05;
elseif CASE==8     % not really part of the set
g_hat 		= 50;
p_hat 		= 0.05;
p_hatIP3 	= 0.05;
end


%% Myosin crossbridge model
global C_Hillmann
K2_c        = 0.5 * C_Hillmann;
K3_c        = 0.4 * C_Hillmann;
K4_c        = 0.1 * C_Hillmann;
K5_c        = 0.5 * C_Hillmann;
K7_c        = 0.1 * C_Hillmann;
gam_cross   = 17 * C_Hillmann;

%% Koningsberger

% nu_r        = 100*inv(0.0075);
% sigp0_r     = 0.0191*inv(0.0075);
% kp_r        = 0.15;
% r0_r        = 20;
% siga0_r     = 1.8e5;
% ka_r        = 0.0006;
% ra_r        = 12;
% rb_r        = 15;
% hb_r        = 3;
% P_r         = 30*inv(0.0075);


%% Kelvin Voigt

P_r         = 4000;  % Pa
rb_r        = 20e-6; % m
h0_r        = 3e-6;  % m 
R0pas_r     = 20e-6;
R0act_r     = 12e-6;
Eact_r      = 233e3;
Epas_r      = 66e3;
nu_r        = 1e4;
