%Run stats simulations with real EEG noise trials and effects
%
%Author: Eric Fields
%Version Date: 22 March 2019

clearvars; close all;

main_dir = MUSim_main_dir();

n_exp  = 2.5e2;
n_perm = 2.5e2;
save_results = false;

effect = fullfile(main_dir, 'data', 'NonCon_N400_restricted.mat');
time_wind = [200, 600];
electrodes = 1:32; %[10, 14, 21, 22, 23];
mult_comp_method = 'clust';

run_real_erp_sim(effect, time_wind, electrodes, mult_comp_method, n_exp, n_perm, save_results);
