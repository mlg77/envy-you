function [ desired_save ] = nvu_script_easy_run_myfun( my_J_PLC )

%% Demonstration script for new NVU model reverting back to version 1.0
%      Author: Michelle Goodman
%      Date: 3/2/2015

%      For reference 
%           "Old" refers to version 1.0 orgional prior to Ca and speedfix 
%           "New" refers to version 1.1 with Ca in astrocyte and speedfix edits

% First copy from old code running types
odeopts = odeset('RelTol', 1e-03, 'AbsTol', 1e-03, 'MaxStep', 1); 
% odeopts = odeset('Vectorized', 1); This was the latest

nv = NVU(Astrocyte(), ...
    WallMechanics(), ...
    SMCEC('J_PLC', my_J_PLC), ...
    'odeopts', odeopts);


%% Changing parameters, initial conditions, simulation time to be the same as the old
nv.smcec.params.J_PLC = my_J_PLC; % new is 0.4;  % (muM s-1) EC agonist concentration  
nv.T = linspace(0, 500, 1000);
nv.simulate();

%desired_save  = [nv.T, 1e6 * nv.out('R'), ];
desired_save  = [nv.T,1e6 * nv.out('R')];

end