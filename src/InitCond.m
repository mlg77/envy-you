function STATES = InitCond()
    % Below the initial conditions of the differential equation are given.
    % They are chosen, such that the system is in steady state at t=0
    
    all_indices()
    
    STATES(ind.R_k)     = 0.061e-6;     %'wi in component wi (metre)'
    STATES(ind.N_Na_k)  = 0.99796e-3;   %'N_Nai in component N_Nai (micromolar_metre)'
    STATES(ind.N_K_k)   = 5.52782e-3;   %'N_Ki in component N_Ki (micromolar_metre)'
    STATES(ind.N_HCO3_k)= 0.58804e-3;   %'N_HCO3i in component N_HCO3i (micromolar_metre)'
    STATES(ind.N_Cl_k)  = 0.32879e-3;   %'N_Cli in component N_Cli (micromolar_metre)'
    STATES(ind.N_Na_s)  = 4.301041e-3;  %'N_Nao in component N_Nao (micromolar_metre)'
    STATES(ind.N_K_s)   = 0.0807e-3;    %'N_Ko in component N_Ko (micromolar_metre)'
    STATES(ind.N_HCO3_s)= 0.432552e-3;  %'N_HCO3o in component N_HCO3o (micromolar_metre)'
    STATES(ind.K_p)     = 3e3;         % uM,  [K+] in de perivascular space
    STATES(ind.w_k)     = 0.1815e-3;    % [-]  BK-Channel open probability
   
    STATES(ind.Ca_i)    = 0.1;            % calcium concentration in cytosol
    STATES(ind.s_i)     = 0.1;            % calcium concentration in sacroplasmatic reticulum
    STATES(ind.v_i)     = -60;            % mV celmembrane of SMC
    STATES(ind.w_i)     = 0.1;            % open state probability of calcium-activated K channels
    STATES(ind.I_i)     = 0.1;            % IP3 concentration
    
    STATES(ind.K_i)     = 100e3;            %uM [K+] in SMC
    
    STATES(ind.Ca_j)    = 0.1;            % calcium concentration in EC cytosol
    STATES(ind.s_j)     = 0.1;            % calcium concentration in endoplasmatic reticulum
    STATES(ind.v_j)     = -75;            % mV celmembrane of EC
    STATES(ind.I_j)     = 0.1;            % IP3 concentration in EC
    
    STATES(ind.Mp)      = 0.25;
    STATES(ind.AMp)     = 0.25;
    STATES(ind.AM)      = 0.25;
    
    STATES(ind.R)       = 15e-6;
    STATES(ind.N_O2)    = 2e-2;
   
    
    
    STATES(ind.N_sa_Em)             = -70   ;    %mV membrane potential of soma
    STATES(ind.N_d_Em)              = -70   ;    %mV membrane potential of dendrite
    STATES(ind.N_sa_Na_NaP_m1_GHK)  = 0.012869;
    STATES(ind.N_sa_Na_NaP_h1_GHK)  = 0.9718;
    STATES(ind.N_sa_K_KDR_m2_GHK)   = 0.0012175;
    STATES(ind.N_sa_K_KA_m3_GHK)    = 0.1193;
    STATES(ind.N_sa_K_KA_h2_GHK)    = 0.12053;
    STATES(ind.N_d_Na_NaP_m4_GHK)   = 0.012869;
    STATES(ind.N_d_Na_NaP_h3_GHK)   = 0.9718;
    STATES(ind.N_d_Na_NMDA_m5_GHK)  = 0.00087339;
    STATES(ind.N_d_Na_NMDA_h4_GHK)  = 0.9898;
    STATES(ind.N_d_K_KDR_m6_GHK)    = 0.0012175;
    STATES(ind.N_d_K_KA_m7_GHK)     = 0.1193;
    STATES(ind.N_d_K_KA_h5_GHK)     = 0.12053;
%     STATES(ind.N_d_K_NMDA_m8_GHK)   = 0.0008;
%     STATES(ind.N_d_K_NMDA_h6_GHK)   = 0.9898;
    STATES(ind.N_sa_Na)             = 10    ;    %mM concentration of Na inside the soma
    STATES(ind.N_sa_K)              = 133.5 ;    %mM concentration of K inside the soma
    STATES(ind.N_sa_Cl)             = 10.4287  ;    %mM concentration of Cl inside the soma
    STATES(ind.N_d_Na)              = 10    ;    %mM concentration of Na inside the dendrite
    STATES(ind.N_d_K)               = 133.5 ;    %mM concentration of K inside the dendrite
    STATES(ind.N_d_Cl)              = 10.4287  ;    %mM concentration of Cl inside the dendrite
    STATES(ind.N_e_K_buffer)        = 134.9701;    %mM free buffer concentration   
    STATES(ind.N_e_Na)              = 140   ;    %mM concentration of Na in the extracellular space    
    STATES(ind.N_e_K)               = 3.5   ;    %mM concentration of K in the extracellular space
    STATES(ind.N_e_Cl)              = 143.5 ;    %mM concentration of Cl in the extracellular space
   
end