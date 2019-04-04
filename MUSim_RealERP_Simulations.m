%Run stats simulations with real EEG noise trials and effects
%
%Author: Eric Fields
%Version Date: 2 April 2019

clearvars; close all;

main_dir = MUSim_main_dir();

%Add EEGLAB, MUT, and FMUT to path if on cluster
if isunix()
    addpath('/gsfs0/data/fields/Documents/MATLAB/eeglab14_1_2b_ECF');
end

n_exp  = 1e2;
n_perm = 1e3;
alpha = 0.05;
save_results = true;

effect = fullfile(main_dir, 'data', 'NonCon_N400_restricted.mat');
factor_levels = 2;
time_wind = [0, 1000];
electrodes = 1:32; %[10, 14, 21, 22, 23];
mult_comp_method = 'clust05';

diary(fullfile(main_dir, 'results', 'MUSim_sirius_test.txt'));
run_real_erp_sim(effect, time_wind, electrodes, factor_levels, n_exp, n_perm, alpha, save_results);
diary off;
