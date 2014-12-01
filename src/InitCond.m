function STATES = InitCond()
    % Below the initial conditions of the differential equation are given.
    % They are chosen, such that the system is in steady state at t=0
    
    all_indices()
    global nonDim

    if strcmp(nonDim,'ON') == 1
        STATES = zeros(size(fields(ind)));
        for i = 1:size(fields(ind)) 
            STATES(i) = 1;
        end
        
    elseif strcmp(nonDim,'OFF') == 1
        STATES(ind.R_k_star)     = 0.061e-6;     %'wi in component wi (metre)'
        STATES(ind.N_Na_k_star)  = 0.99796e-3;   %'N_Nai in component N_Nai (micromolar_metre)'
        STATES(ind.N_K_k_star)   = 5.52782e-3;   %'N_Ki in component N_Ki (micromolar_metre)'
        STATES(ind.N_HCO3_k_star)= 0.58804e-3;   %'N_HCO3i in component N_HCO3i (micromolar_metre)'
        STATES(ind.N_Cl_k_star)  = 0.32879e-3;   %'N_Cli in component N_Cli (micromolar_metre)'
        STATES(ind.N_Na_s_star)  = 4.301041e-3;  %'N_Nao in component N_Nao (micromolar_metre)'
        STATES(ind.N_K_s_star)   = 0.0807e-3;    %'N_Ko in component N_Ko (micromolar_metre)'
        STATES(ind.N_HCO3_s_star)= 0.432552e-3;  %'N_HCO3o in component N_HCO3o (micromolar_metre)'
        STATES(ind.K_p_star)     = 3e3;         % uM,  [K+] in the perivascular space
        STATES(ind.w_k_star)     = 0.1815e-3;    % [-]  BK-Channel open probability



        STATES(ind.Ca_i_star)    = 0.1;            % calcium concentration in cytosol
        STATES(ind.s_i_star)     = 0.1;            % calcium concentration in sacroplasmatic reticulum
        STATES(ind.v_i_star)     = -60;            % mV celmembrane of SMC
        STATES(ind.w_i_star)     = 0.1;            % open state probability of calcium-activated K channels
        STATES(ind.I_i_star)     = 0.1;            % IP3 concentration

        STATES(ind.K_i_star)     = 100e3;            %uM [K+] in SMC

        STATES(ind.Ca_j_star)    = 0.1;            % calcium concentration in EC cytosol
        STATES(ind.s_j_star)     = 0.1;            % calcium concentration in endoplasmatic reticulum
        STATES(ind.v_j_star)     = -75;            % mV celmembrane of EC
        STATES(ind.I_j_star)     = 0.1;            % IP3 concentration in EC

        STATES(ind.Mp_star)      = 0.25;
        STATES(ind.AMp_star)     = 0.25;
        STATES(ind.AM_star)      = 0.25;

        STATES(ind.R_star)       = 15e-6;

%               %Hannah:
%         STATES(ind.Ca_k)      =0.05e-3;       %uM Bennet 2008
%         STATES(ind.s_k)      =0.1e-3;
%         STATES(ind.h_k)      =0.1e-3;
%         STATES(ind.I_k)      =0.01e-3;       %uM Bennet 2008 0.16 uM according to Hadvield David
%         STATES(ind.eet_k)    =0.1e-3;
    end
        
end