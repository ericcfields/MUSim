%Run stats simulations with real EEG noise trials and effects
%
%Author: Eric Fields
%Version Date: 5 April 2019

%% SET-UP

clearvars; close all;

%Get main project directory
main_dir = MUSim_main_dir();

%Add EEGLAB, MUT, and FMUT to path if on cluster
if isunix()
    addpath('/gsfs0/data/fields/Documents/MATLAB/eeglab14_1_2b_ECF');
end


%% SIMULATION PARAMETERS

%noise trials file
noise = fullfile(main_dir, 'data', 'noise_trials.mat');

%Simulated data parameters
n_exp  = 1e3;
n_perm = 1e3;
n_trials = 20;
n_subs = 24;
cond_trials = 20;
error_mult = 1;
ind_var_factor = 0.1;

%Effect to simulate (.mat file with effect or 'null')
effect = fullfile(main_dir, 'data', 'NonCon_N400_reduced.mat');
factor_levels = 2;

%Analysis parameters
time_wind = [300 500];
electrodes = [10, 14, 21, 22, 23];
alpha = 0.05;

%File for saving results
output_file = false; %sfullfile(main_dir, 'results', 'MUSim_results.txt');


%% RUN SIMULATION

run_real_erp_sim(noise, effect, time_wind, electrodes, factor_levels, n_exp, n_perm, n_subs, cond_trials, error_mult, ind_var_factor, alpha, output_file);

