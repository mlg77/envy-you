function [dy] = DEsyst (time,state)
% Below all conservation equation are calculated.
% Note that "getRef" is the function that contains the input stimulus of
% the model.
global Ca_0 v_0 K_0 R_0 Rk_0 T_c

dy = zeros(size(state));
all_constants(); % All constants used in this model
all_indices() ; % All indices used in this model

% All additional equation are calculated using the state variables. They are stored in three matrices; AC, SMC and EC
% Note that EC is empty at the moment.

[NE,AC,SMC,EC] = all_fluxes(time, state);

% Astrocyte
dy(ind.R_k_star) = T_c*( L_p * (AC(flu.Na_k) + AC(flu.K_k) + AC(flu.Cl_k) + AC(flu.HCO3_k) - AC(flu.Na_s) - AC(flu.K_s) - AC(flu.Cl_s) - AC(flu.HCO3_s) + X_k / (state(ind.R_k_star)*Rk_0))); % m s-1
dy(ind.N_Na_k_star) = T_c*( -AC(flu.J_Na_k) - 3 * AC(flu.J_NaK_k) + AC(flu.J_NKCC1_k) + AC(flu.J_NBC_k )); % uMm s-1
dy(ind.N_K_k_star) = T_c*( -AC(flu.J_K_k ) + 2 * AC(flu.J_NaK_k) + AC(flu.J_NKCC1_k) + AC(flu.J_KCC1_k) - AC(flu.J_BK_k)); % uMm s-1
dy(ind.N_HCO3_k_star) = T_c*( 2 * AC(flu.J_NBC_k)); % uMm s-1
dy(ind.N_Cl_k_star) = T_c*( dy(ind.N_Na_k_star) + dy(ind.N_K_k_star) - dy(ind.N_HCO3_k_star)); % uMm s-1, modified equation compared to the one of Ostby
% dy(ind.N_Na_s_star) = T_c*( -k_C * getRef(time,'ft') - dy(ind.N_Na_k_star)); % uMm s-1
dy(ind.N_Na_s_star) = T_c*( -k_C * 0 - dy(ind.N_Na_k_star)); % uMm s-1



% dy(ind.N_K_s_star) = T_c*(    k_C * getRef(time,'ft') - dy(ind.N_K_k_star) +AC(flu.J_BK_k)) ; % uMm s-1
dy(ind.N_K_s_star) = T_c*(    k_C * 0 - dy(ind.N_K_k_star) +AC(flu.J_BK_k)) ; % uMm s-1


dy(ind.N_HCO3_s_star) = T_c*(  - dy(ind.N_HCO3_k_star)); % uMm s-1
dy(ind.K_p_star) = T_c*(      AC(flu.J_BK_k) / (VR_pa*(state(ind.R_k_star)*Rk_0)) + (SMC(flu.J_KIR_i))/(VR_ps)); % uM s-1
dy(ind.w_k_star) = T_c*(      AC(flu.phi_w) * (AC(flu.w_inf) - state(ind.w_k_star))); % s-1

% %astr. DE for Ca2+ (Hannah)
%dy(ind.ck)=T_c*( B_cyt*(AC(flu.J_ip3)-AC(flu.J_pump)+AC(flu.J_ERleak))); % FARR astrocytic calcium concentration
dy(ind.Ca_k_star)= T_c*(         AC(flu.B_cyt)*(AC(flu.J_ip3)-AC(flu.J_pump)+AC(flu.J_ERleak))); %L&E astrocytic calcium concentration
dy(ind.s_k_star)= T_c*(         -1/VR_ERcyt*dy(ind.Ca_k_star)); % ER calcium concentration   !!!!!!!!!!!!!!!!!!!!!
dy(ind.h_k_star)= T_c*(         k_on*(k_inh-(AC(flu.ck)+k_inh)*AC(flu.hk))); %the action of the IP3 receptors that have not been inactivated by Calcium
dy(ind.I_k_star)= T_c*(         r_h*AC(flu.G_pr)-k_deg*AC(flu.ik)); %IP3 concentration
%dy(ind.eetk)=V_eet*(AC(flu.ck)-ck_min)-k_eet*AC(flu.eetk); % FARR EET concentration

if AC(flu.ck)> ck_min
dy(ind.eet_k_star)= T_c*( V_eet*(AC(flu.ck)-ck_min)-k_eet*AC(flu.eetk)); % L&E EET concentration
else
dy(ind.eet_k_star)= T_c*( -k_eet*AC(flu.eetk)); %L&E
end

% Smooth muscle cell
dy(ind.Ca_i_star) = T_c*( SMC(flu.Ca_coup_i) + SMC(flu.rho_i) * (SMC(flu.J_CICR_i) + SMC(flu.J_IP3_i) + SMC(flu.J_leak_i) - SMC(flu.J_SRuptake_i) - SMC(flu.J_extrusion_i) - SMC(flu.J_VOCC_i) + SMC(flu.J_NaCa_i) + 0.1*SMC(flu.J_stretch_i)));
dy(ind.s_i_star) = T_c*( -SMC(flu.J_CICR_i) - SMC(flu.J_leak_i) + SMC(flu.J_SRuptake_i) );
dy(ind.v_i_star) = T_c*( SMC(flu.v_coup_i) + gam * (-SMC(flu.J_NaK_i) - SMC(flu.J_Cl_i) - 2*SMC(flu.J_VOCC_i) - SMC(flu.J_NaCa_i) - SMC(flu.J_K_i) - SMC(flu.J_stretch_i) - SMC(flu.J_KIR_i)));
dy(ind.w_i_star) = T_c*( lab * (SMC(flu.Kactivation_i) - state(ind.w_i_star)));
dy(ind.I_i_star) = T_c*( SMC(flu.IP3_coup_i) - SMC(flu.J_degrad_i)) ;
dy(ind.K_i_star) = T_c*( -SMC(flu.J_KIR_i) - SMC(flu.J_K_i) + SMC(flu.J_NaK_i)); % uM s-1

% Endothelium cell
dy(ind.Ca_j_star) = T_c*( EC(flu.Ca_coup_j) + EC(flu.rho_j) * (EC(flu.J_IP3_j) - EC(flu.J_ERuptake_j) + EC(flu.J_CICR_j) - EC(flu.J_extrusion_j) + EC(flu.J_leak_j) + EC(flu.J_cation_j) + EC(flu.J_0_j) + EC(flu.J_stretch_j)));
dy(ind.s_j_star) = T_c*( EC(flu.J_ERuptake_j) - EC(flu.J_CICR_j) - EC(flu.J_leak_j) );
dy(ind.v_j_star) = T_c*( - 1/C_m * ( EC(flu.J_K_j)	+ EC(flu.J_R_j) ) + EC(flu.v_coup_j));
dy(ind.I_j_star) = T_c*( EC(flu.IP3_coup_j) + J_PLC - EC(flu.J_degrad_j)) ;

% Myosin crossbridge model
dy(ind.Mp_star) = T_c*( K4_c*state(ind.AMp_star) + SMC(flu.K1_c) *SMC(flu.M) - (K2_c + K3_c)*state(ind.Mp_star));
dy(ind.AMp_star) = T_c*( K3_c*state(ind.Mp_star) + SMC(flu.K6_c) *state(ind.AM_star) - (K4_c + K5_c)*state(ind.AMp_star)); % K7_c was corrected to K4_c
dy(ind.AM_star) = T_c*( K5_c*state(ind.AMp_star) - (K7_c + SMC(flu.K6_c) )*state(ind.AM_star));

% Radius change
%Kelvin Voigt
F_r =               state(ind.AMp_star) + state(ind.AM_star);
E_r =               Epas_r + F_r*(Eact_r -Epas_r);
R0_r=               R0pas_r + F_r*R0pas_r*(0.6 - 1);
dy(ind.R_star)= T_c*( R0pas_r/nu_r *((state(ind.R_star)*R_0)*P_r/SMC(flu.h_r) - E_r * (((state(ind.R_star)*R_0) - R0_r)/R0_r)));


% if F_r1 <= 0.4
% F_r = 0.4;
% else
% F_r = F_r1;
% end


%dy(ind.R) = 1/nu_r *( R0pas_r*(state(ind.R_star)*R_0)*P_r /SMC(flu.h_r) - Epas_r * ( (state(ind.R_star)*R_0) - R0pas_r)...
% - F_r/0.8 * ((Etot_r - Epas_r) * (state(ind.R_star)*R_0) + Epas_r*R0pas_r - Etot_r*R0act_r));


%Koningsberger
% siga_r = siga0_r *((state(AMp) + state(AM))/0.8) * exp(-ka_r *(state(R)+h_r-ra_r)^2);
%

% if state(R) >= r0_r
% sigp_r = sigp0_r * (exp(kp_r * (state(R)-r0_r)) -1);
% else
% sigp_r = sigp0_r * kp_r * (1 - (state(R)^2/r0_r^2)^(-3/2));
% end
%
%
% dy(R) = 1/nu_r * ( P_r*state(R)/h_r - sigp_r - siga_r);
end