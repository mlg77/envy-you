% New script - plot some graphs as a check to see everything is working
% with new model

clear;
odeopts = odeset('Vectorized', 1);

% New model with neuron class:
nv = NVU2(Astrocyte2(), WallMechanics(), SMCEC(), Neuron(), ...
    'odeopts', odeopts);

% Old model without neuron class:
%nv = NVU(Astrocyte(), WallMechanics(), SMCEC(), ...
    %'odeopts', odeopts);

nv.simulate();

figure(1);
plot(nv.T, nv.out('N_Na_n'));
title('Na in neuron'); ylim([-1e-5,2e-3]); xlim([0,500]);
figure(2);
plot(nv.T, nv.out('N_K_n'));
title('K in neuron'); ylim([-1e-5,2e-3]); xlim([0,500]);
figure(3);
plot(nv.T, nv.out('N_Na_s'));
title('Na in SC');
figure(4);
plot(nv.T, nv.out('N_K_s'));
title('K in SC');
figure(5);
plot(nv.T, nv.out('R'));
figure(6);
plot(nv.T, nv.out('J_Na_n'));