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

n_exp  = 1e3;
n_perm = 1e3;
n_trials = 20;
alpha = 0.05;
save_results = true;

effect = 'null'; %fullfile(main_dir, 'data', 'namesP300_reduced.mat');
factor_levels = 3;
time_wind = [500 750];
electrodes = [10, 14, 21, 22, 23];

run_real_erp_sim(effect, time_wind, electrodes, factor_levels, n_exp, n_perm, n_trials, alpha, save_results);
