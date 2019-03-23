%Run stats simulations with real EEG noise trials and effects
%
%Author: Eric Fields
%Version Date: 22 March 2019

clearvars; close all;

main_dir = MUSim_main_dir();

%Add EEGLAB, MUT, and FMUT to path if on cluster
if isunix()
    addpath('/gsfs0/data/fields/Documents/MATLAB/eeglab14_1_1b_ECF'); %NEED TO CHANGE
    addpath('/gsfs0/data/fields/Documents/MATLAB/FMUT_0.3.5');
    addpath('/gsfs0/data/fields/Documents/MATLAB/dmgroppe-Mass_Univariate_ERP_Toolbox-10dc5c7');
end

n_exp  = 5e3;
n_perm = 5e3;
save_results = false;

effect = fullfile(main_dir, 'data', 'NonCon_N400_restricted.mat');
factor_levels = 2;
time_wind = [300, 500];
electrodes = [10, 14, 21, 22, 23];
mult_comp_method = 'Fmax';

diary(fullfile(main_dir, 'results', 'MUSim_sirius_test.txt'));
tic;
run_real_erp_sim(effect, time_wind, electrodes, factor_levels, mult_comp_method, n_exp, n_perm, save_results);
toc;
diary off;
