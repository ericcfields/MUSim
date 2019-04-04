%Function to test Type I and Type II error rates for mass univariate approaches
%with data constructed from real EEG noise trials and real ERP effects
%
%Author: Eric Fields
%Version Date: 4 April 2019

function run_real_erp_sim(effect, time_wind, electrodes, factor_levels, n_exp, n_perm, alpha, save_results)

    
    %% ####################################################################
    %####################### PARAMETERS AND SET-UP ########################
    %######################################################################
    
    main_dir = MUSim_main_dir();
    
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; %#ok<ASGLU>
    close all;
    rmpath('C:\Users\ecfne\Documents\MATLAB\eeglab14_1_2b_ECF\plugins\FMUT_0.4.1');
    addpath('C:\Users\ecfne\Documents\Eric\Coding\FMUT_development\FMUT');
    
    %% ~~~~~ SIMULATION PARAMETERS ~~~~~
    
    global VERBLEVEL
    VERBLEVEL = 0;
    
    %Load template GND for time point and electrode information
    load(fullfile(main_dir, 'data', 'MUSim_template.GND'), '-mat');
    %Load EEG noise trials
    load(fullfile(main_dir, 'data', 'AXCPT_noise_trials_128Hz_30Hz_49subs.mat'));
    %Load effect data
    if effect
        load(effect);
        if ~isequal(factor_levels, size(effects_data, 3)) %#ok<NODEF>
            error('factor_levels input doesn''t match effects data');
        end
    end
    
    %Simulated data characteristics
    n_subs = 24;
    n_time_pts = length(GND.time_pts);
    bins = 1:prod(factor_levels);
    error_mult  = ones(length(bins), 1) * 3;
    cond_trials = ones(length(bins), 1) * 20;
    sub_param = [1 0.1];
    
    %Parameters for permutation ANOVA
    chan_hood = 75;
    chan_hood = spatial_neighbors(GND.chanlocs(electrodes), chan_hood, []);
    thresh_p = 0.05;
    
    %Parameters for parametric ANOVA
    param_time_wind  = time_wind;
    param_electrodes = electrodes;
    [~, param_start_sample] = min(abs(GND.time_pts - param_time_wind(1)));
    [~, param_end_sample  ] = min(abs(GND.time_pts - param_time_wind(2)));
    
    %Other useful numbers
    n_electrodes = length(electrodes);
    n_conds      = length(bins);
    n_trials     = sum(cond_trials);
    
 
    %% ####################################################################
    %#################### SIMULTATIONS ####################################
    %######################################################################
    
    %% ~~~~~ SET-UP ~~~~~
    
    %Find relevant time points for time window
    [~, start_sample] = min(abs(GND.time_pts - time_wind(1)));
    [~, end_sample  ] = min(abs(GND.time_pts - time_wind(2)));
    n_sample_time_pts = end_sample - start_sample + 1;
    
    %Pre-allocate variables
    h_uncorr  = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_Fmax    = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_clust05 = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_clust01 = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_bh      = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_by      = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_bky     = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_mean_amp = NaN(n_exp, 1);
    p_uncorr  = NaN(n_exp, n_electrodes, n_sample_time_pts);
    p_Fmax    = NaN(n_exp, n_electrodes, n_sample_time_pts);
    p_clust05 = NaN(n_exp, n_electrodes, n_sample_time_pts);
    p_clust01 = NaN(n_exp, n_electrodes, n_sample_time_pts);
    p_bh      = NaN(n_exp, n_electrodes, n_sample_time_pts);
    p_by      = NaN(n_exp, n_electrodes, n_sample_time_pts);
    p_mean_amp = NaN(n_exp, 1);

    %Randomly select subjects to use in each experiment below
    sub_sample = NaN(n_exp, n_subs);
    for i = 1:n_exp
        sub_sample(i, :) = randsample(1:size(all_noise_trials, 1), n_subs); %#ok<USENS>
    end
    
    %Load effects data
    if effect
        effects_data = effects_data(electrodes,:,:);
    else
        effects_data = [];
    end
    

    %% ~~~~~ SIMULATE EXPERIMENTS ~~~~~
    tic
    %Conduct n_exp simulated experiments
    parfor i = 1:n_exp
        
        %% ~~~~~ GENERATE SIMULATED DATA ~~~~~
        
        %Randomly select trials from each subject
        sim_data = NaN(n_electrodes, n_time_pts, n_conds, n_subs);
        for s = 1:n_subs
            sim_data(:,:,:,s) = mean(reshape(all_noise_trials{sub_sample(i,s),2}(electrodes, :, randsample(1:size(all_noise_trials{sub_sample(i,s),2},3), n_trials)), n_electrodes, n_time_pts, n_conds, []), 4);
        end
        
        %Add effects if any
        if effect
            %Generate subject effects for this experiment
            sub_effects = normrnd(sub_param(1), sub_param(2), n_subs, n_conds);
            for s = 1:n_subs
                for b = 1:n_conds
                    sim_data(:, :, b, s) = sim_data(:, :, b, s)*error_mult(b) + effects_data(:, :, b)*sub_effects(s, b);
                end
            end
        end
        
        %% ~~~~~ RUN STATS ~~~~~
        
        %Get analysis data
        data = sim_data(:, start_sample:end_sample, bins, :);
        data = reshape(data,[n_electrodes, n_sample_time_pts, factor_levels, n_subs]);
        
        %Run ANOVA
        [F_obs, F_dist, df_effect, df_res] = perm_rbANOVA(data, 3, n_perm);
        
        %Fmax
        [h_Fmax(i, :, :), p_Fmax(i, :, :)] = Fmax_corr(F_obs, F_dist, alpha);
        thresh_F = finv(1-0.05, df_effect, df_res);
        
        %cluster
        [h_clust05(i, :, :), p_clust05(i, :, :)] = Fclust_corr(F_obs, F_dist, alpha, chan_hood, thresh_F);
        thresh_F = finv(1-0.01, df_effect, df_res);
        [h_clust01(i, :, :), p_clust01(i, :, :)] = Fclust_corr(F_obs, F_dist, alpha, chan_hood, thresh_F);
        
        %FDR
        p_uncorr(i, :, :) = 1 - fcdf(F_obs, df_effect, df_res);
        h_uncorr(i, :, :) = p_uncorr(i, :, :) <= alpha;
        [h_bh(i, :, :), ~, ~, p_bh(i, :, :)] = fdr_bh(p_uncorr(i, :, :), alpha, 'pdep', 'no');
        [h_by(i, :, :), ~, ~, p_by(i, :, :)] = fdr_bh(p_uncorr(i, :, :), alpha, 'dep', 'no');
        h_bky(i, :, :) = fdr_bky(p_uncorr(i, :, :), alpha, 'no');
        
        %Mean amplitude ANOVA
        if isequal(electrodes, 1:32)
            mean_wind_data = mean(mean(sim_data(param_electrodes, param_start_sample:param_end_sample, bins, :), 1), 2);
            mean_wind_data = reshape(mean_wind_data, [size(mean_wind_data, 1), size(mean_wind_data, 2), factor_levels, size(mean_wind_data, 4)]);
        elseif isequal(electrodes, param_electrodes)
            mean_wind_data = mean(mean(sim_data(:, param_start_sample:param_end_sample, bins, :), 1), 2);
            mean_wind_data = reshape(mean_wind_data, [size(mean_wind_data, 1), size(mean_wind_data, 2), factor_levels, size(mean_wind_data, 4)]);
        else
            error('Can''t determine electrodes for parametric test');
        end
        param_results = calc_param_ANOVA(mean_wind_data, [], (1:length(factor_levels > 1))+2, alpha);
        p_mean_amp(i) = param_results.p;
        h_mean_amp(i) = param_results.h;
        
    end
    toc
    
    %Calcluate Dunn-Sidak corrected results
    p_sidak = 1-(1-p_uncorr).^(n_electrodes*n_time_pts);
    p_sidak(p_sidak>1) = 1;
    h_sidak = p_sidak < alpha;
    
    
    %% ####################################################################
    %################# SAVE AND REPORT RESULTS ############################
    %######################################################################
    
    %Save results
    time_stamp = sprintf('%s_%s', datestr(datetime('now'),'ddmmmyy'), datestr(datetime('now'),'HHMM'));
    if save_results
        %save(sprintf('R:\\Public\\GK_lab\\Eric\\Stats Simulations\\erps\\RealERP_simulation_%s.mat', time_stamp), '-v7.3'); 
        diary(fullfile(main_dir, 'results', 'MUSim_results.txt'));
    end
    
    if ~effect
        effect = 'none';
        effect_description = 'all effects null';
    end
    
    %Print output
    fprintf('\n');
    fprintf('----------------------------------------------------------------------------------\n')
    fprintf('****SIMULATION SUMMARY****\n')
    fprintf('%s\n\n', time_stamp);
    
    fprintf('Using real EEG noise trials\n');
    fprintf('%d simulated experiments, %d permutations each, %d subjects\n\n', n_exp, n_perm, n_subs);
    fprintf('\nEffect: %s', effect);
    fprintf('\nEffect description: %s\n', effect_description);
    fprintf('Time window: %d - %d\n', time_wind(1), time_wind(2));
    fprintf('Electrodes: ');
    fprintf([sprintf('%d, ', electrodes(1:end-1)), num2str(electrodes(end))]);
    fprintf('\nError multiplier = ');
    fprintf('%.1f  ', error_mult);
    fprintf('\nFactor levels: ');
    fprintf('%d  ', factor_levels);
    fprintf('\nTrials = ');
    fprintf('%d  ', cond_trials);
    fprintf('\n');
    
    fprintf('\nCLUSTER MASS 0.05 RESULTS\n')
    fprintf('Rejection rate all time points (individually) = %f\n', mean(h_clust05(:)));
    fprintf('Family-wise rejection rate = %f\n', mean(any(reshape(h_clust05, n_exp, n_electrodes*n_sample_time_pts)')));
    
    fprintf('\nCLUSTER MASS 0.01 RESULTS\n')
    fprintf('Rejection rate all time points (individually) = %f\n', mean(h_clust01(:)));
    fprintf('Family-wise rejection rate = %f\n', mean(any(reshape(h_clust01, n_exp, n_electrodes*n_sample_time_pts)')));
    
    fprintf('\nFMAX RESULTS\n')
    fprintf('Rejection rate all time points (individually) = %f\n', mean(h_Fmax(:)));
    fprintf('Family-wise rejection rate = %f\n', mean(any(reshape(h_Fmax, n_exp, n_electrodes*n_sample_time_pts)')));
    
    fprintf('\nBH FDR RESULTS\n')
    fprintf('Rejection rate all time points (individually) = %f\n', mean(h_bh(:)));
    fprintf('Family-wise rejection rate = %f\n', mean(any(reshape(h_bh, n_exp, n_electrodes*n_sample_time_pts)')));
    
    fprintf('\nBY FDR RESULTS\n')
    fprintf('Rejection rate all time points (individually) = %f\n', mean(h_by(:)));
    fprintf('Family-wise rejection rate = %f\n', mean(any(reshape(h_by, n_exp, n_electrodes*n_sample_time_pts)')));
    
    fprintf('\nBKY FDR RESULTS\n')
    fprintf('Rejection rate all time points (individually) = %f\n', mean(h_bky(:)));
    fprintf('Family-wise rejection rate = %f\n', mean(any(reshape(h_bky, n_exp, n_electrodes*n_sample_time_pts)')));
    
    fprintf('\nMEAN WINDOW/REGION PARAMETRIC F-TEST RESULTS\n')
    fprintf('Rejection rate = %f\n\n', mean(h_mean_amp));
    
    fprintf('DUNN-SIDAK CORRECTED PARAMETRIC F-TEST RESULTS\n');
    fprintf('Rejection rate all time points (individually) = %f\n', mean(h_sidak(:)));
    fprintf('Family-wise rejection rate = %f\n', mean(any(reshape(h_sidak, n_exp, n_electrodes*n_sample_time_pts)')));
    
    fprintf('\nUNCORRECTED PARAMETRIC F-TEST RESULTS\n');
    fprintf('Rejection rate all time points (individually) = %f\n', mean(h_uncorr(:)));
    fprintf('Family-wise rejection error rate = %f\n', mean(any(reshape(h_uncorr, n_exp, n_electrodes*n_sample_time_pts)')));
    
    fprintf('\n----------------------------------------------------------------------------------\n\n')
    
    diary off
    
end
