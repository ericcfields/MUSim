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

n_exp  = 1e1;
n_perm = 1e2;
save_results = false;

effect = fullfile(main_dir, 'data', 'NonCon_N400_restricted.mat');
factor_levels = 2;
time_wind = [0, 1000];
electrodes = 1:32; %[10, 14, 21, 22, 23];
mult_comp_method = 'clust05';

diary(fullfile(main_dir, 'results', 'MUSim_sirius_test.txt'));
tic;
run_real_erp_sim(effect, time_wind, electrodes, factor_levels, mult_comp_method, n_exp, n_perm, save_results);
toc;
diary off;
