%Run stats simulations with real EEG noise trials and effects
%
%Author: Eric Fields
%Version Date: 2 July 2019
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
n_exp  = 1e2; %number of simulated experiments
n_perm = 5e2; %permutations per experiment for Fmax and clust procedures
error_mult = 1;   %factor to multiple error standard deviation by (can be array for testing unequal variances)
ind_var_factor = 0.1; %standard deviation of multiplier for individual differences in effects

%Analysis parameters
alpha = 0.05;

%File for saving results
text_output = fullfile(main_dir, 'results', 'MUSim_null_results.txt');


%% RUN SIMULATIONS

effect = 'null';
factor_levels = [3, 3];
dims = [3, 4];
time_windows = {[0 300], [300 1000]};
electrodes = 1:32;

for n_subs = [40, 25, 16, 12, 8]
    for cond_trials = [40, 20, 10]
        for t = 1:length(time_windows)
            time_wind = time_windows{t};
            mat_output = fullfile(main_dir, 'results', sprintf('MUSim_null_%d-%d_simulation_results.mat', time_wind(1), time_wind(2)));
            run_real_erp_sim(noise, effect, time_wind, electrodes, factor_levels, 3, n_exp, n_perm, n_subs, cond_trials, error_mult, ind_var_factor, alpha, text_output, mat_output);
            return;
        end
    end
end
