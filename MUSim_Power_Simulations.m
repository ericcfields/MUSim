%Run stats simulations with real EEG noise trials and effects
%
%Author: Eric Fields
%Version Date: 17 June 2019
%
%Copyright (c) 2019, Eric C. Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.
%This software is provided "as is" and any express or implied warranties are disclaimed. 

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
n_exp  = 1e4; %number of simulated experiments
n_perm = 5e3; %permutations per experiment for Fmax and clust procedures
n_subs = 24;  %number of subjects in each simulated experiment
cond_trials = 20; %number of trials in each condition
error_mult = 1;   %factor to multiple error standard deviation by (can be array for testing unequal variances)
ind_var_factor = 0.1; %standard deviation of multiplier for individual differences in effects

%Analysis parameters
param_file = fullfile(main_dir, 'MUSim_power_sim_params.csv');
sim_list = readtable(param_file);
alpha = 0.05;

%File for saving results
text_output = fullfile(main_dir, 'results', 'MUSim_power_results.txt');


%% RUN SIMULATIONS

for s = 1:size(sim_list, 1)
    
    if ~sim_list{s, 'include'}
        continue;
    end
    
    if ~strcmpi(sim_list{s, 'effect'}{1}, 'null')
        effect = fullfile(main_dir, 'data', sim_list{s, 'effect'}{1});
    else
        effect = sim_list{s, 'effect'}{1};
    end
    factor_levels = sim_list{s, 'factor_levels'};
    time_wind = [sim_list{s, 'start_time'}, sim_list{s, 'end_time'}];
    electrodes = eval(['[' sim_list{s, 'electrodes'}{1} ']']);
    
    [~, effect_name] = fileparts(effect);
    mat_output = fullfile(main_dir, 'results', sprintf('%s_%d-%d_simulation_results.mat', effect_name, time_wind(1), time_wind(2)));
    
    run_real_erp_sim(noise, effect, time_wind, electrodes, factor_levels, 3, n_exp, n_perm, n_subs, cond_trials, error_mult, ind_var_factor, alpha, text_output, mat_output)
    
end
