%% clear
clearvars, close all, clc

%% cylinderwake steadystate
cylinderwake_steadystate();
cylinderwake_steadystate(1, 50, 5);
cylinderwake_steadystate(2, 50, 5);
cylinderwake_steadystate(3, 50, 5);


%% cylinderwake steadystate bccontrol
cylinderwake_steadystate_bccontrol();
cylinderwake_steadystate_bccontrol(1, 50, 5, 1e-4); 
cylinderwake_steadystate_bccontrol(2, 50, 5, 1e-4); 
cylinderwake_steadystate_bccontrol(3, 50, 5, 1e-4); 


%% cylinderwake tdp pout vout
cylinderwake_tdp_pout_vout()
cylinderwake_tdp_pout_vout(1, 50, 5, 0.0, 1.5, 4096);
cylinderwake_tdp_pout_vout(2, 50, 5, 0.0, 1.5, 4096); 
cylinderwake_tdp_pout_vout(3, 50, 5, 0.0, 1.5, 4096); 


%% cylinderwake tdp pout vout bccontrol
cylinderwake_tdp_pout_vout_bccontrol()
cylinderwake_tdp_pout_vout_bccontrol(1, 50, 5, 1e-4, 3.5, 0.0, 1.5, 4096);
cylinderwake_tdp_pout_vout_bccontrol(1, 50, 5, 1e-4, 3.5, 0.0, 1.5, 4096); 
cylinderwake_tdp_pout_vout_bccontrol(1, 50, 5, 1e-4, 3.5, 0.0, 1.5, 4096); 


%% drivencavity steadystate
drivencavity_steadystate()
drivencavity_steadystate(10, 500, 5); 
drivencavity_steadystate(20, 500, 5);
drivencavity_steadystate(30, 500, 5);


%% drivencavity tdp pout vout distrcontrol
drivencavity_tdp_pout_vout_distcontrol() 
drivencavity_tdp_pout_vout_distcontrol(10, 500, 5, 3.5, 0.0, 1.5, 4096); 
drivencavity_tdp_pout_vout_distcontrol(10, 500, 5, 3.5, 0.0, 1.5, 4096); 
drivencavity_tdp_pout_vout_distcontrol(10, 500, 5, 3.5, 0.0, 1.5, 4096); 


%% cylinderwake stokes
stokes_ss()
stokes_ss(1, 50, 5); 
stokes_ss(2, 50, 5): 
stokes_ss(3, 50, 5); 
