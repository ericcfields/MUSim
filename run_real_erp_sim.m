%Function to test Type I and Type II error rates for mass univariate approaches
%with data constructed from real EEG noise trials and real ERP effects
%
%Author: Eric Fields
%Version Date: 4 April 2019

function run_real_erp_sim(effect, time_wind, electrodes, factor_levels, n_exp, n_perm, n_trials, alpha, save_results)

    
    %% ####################################################################
    %####################### PARAMETERS AND SET-UP ########################
    %######################################################################
    
    main_dir = MUSim_main_dir();
    
    %Start EEGLAB
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; %#ok<ASGLU>
    close all;
    
    %% ~~~~~ SIMULATION PARAMETERS ~~~~~
    
    global VERBLEVEL
    VERBLEVEL = 0;
    
    %Load template GND for time point and electrode information
    load(fullfile(main_dir, 'data', 'MUSim_template.GND'), '-mat', 'GND');
    %Load EEG noise trials
    load(fullfile(main_dir, 'data', 'AXCPT_noise_trials_128Hz_30Hz_49subs.mat'), 'all_noise_trials');
    %Load effect data
    if exist(effect, 'file') && ~strcmpi(effect, 'null')
        load(effect, 'effects_data', 'effect_description');
        if ~isequal(factor_levels, size(effects_data, 3))
            error('factor_levels input doesn''t match effects data');
        end
    else
        effect = false;
    end
    
    %Simulated data characteristics
    n_subs = 24;
    n_time_pts = length(GND.time_pts);
    bins = 1:prod(factor_levels);
    error_mult  = ones(length(bins), 1); %Multiplier applied to standard deviation of each bin
    cond_trials = ones(length(bins), 1) * n_trials; %Number of trials in each bin
    ind_var_factor = 0.1; %SD of normal distribuiton for individual difference multiplier
    
    %Channel neighbor information for cluster tests
    chan_hood = 75;
    chan_hood = spatial_neighbors(GND.chanlocs(electrodes), chan_hood, []);
    
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

    %Randomly select subjects to use in each experiment below
    sub_sample = NaN(n_exp, n_subs);
    for i = 1:n_exp
        sub_sample(i, :) = randsample(1:size(all_noise_trials, 1), n_subs);
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
            sub_data = all_noise_trials{sub_sample(i,s),2}; %#ok<PFBNS,PFIIN>
            sim_data(:,:,:,s) = mean(reshape(sub_data(electrodes, :, randsample(1:size(sub_data,3), n_trials)), n_electrodes, n_time_pts, n_conds, []), 4);
        end
        
        %Add effects if any
        if effect
            %Generate subject effects for this experiment
            sub_effects = normrnd(1, ind_var_factor, n_subs, n_conds);
            for s = 1:n_subs
                for b = 1:n_conds
                    sim_data(:, :, b, s) = sim_data(:, :, b, s)*error_mult(b) + effects_data(:, :, b)*sub_effects(s, b); %#ok<PFBNS>
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
        h_Fmax(i, :, :) = Fmax_corr(F_obs, F_dist, alpha);
        thresh_F = finv(1-0.05, df_effect, df_res);
        
        %cluster
        h_clust05(i, :, :) = Fclust_corr(F_obs, F_dist, alpha, chan_hood, thresh_F);
        thresh_F = finv(1-0.01, df_effect, df_res);
        h_clust01(i, :, :) = Fclust_corr(F_obs, F_dist, alpha, chan_hood, thresh_F);
        
        %FDR
        p_uncorr(i, :, :) = 1 - fcdf(F_obs, df_effect, df_res);
        h_uncorr(i, :, :) = p_uncorr(i, :, :) <= alpha;
        h_bh(i, :, :)  = fdr_bh(p_uncorr(i, :, :), alpha, 'pdep', 'no');
        h_by(i, :, :)  = fdr_bh(p_uncorr(i, :, :), alpha, 'dep', 'no');
        h_bky(i, :, :) = fdr_bky(p_uncorr(i, :, :), alpha, 'no');
        
        %Mean amplitude ANOVA
        mean_wind_data = mean(mean(data, 1), 2);
        mean_wind_data = reshape(mean_wind_data, [size(mean_wind_data, 1), size(mean_wind_data, 2), factor_levels, size(mean_wind_data, 4)]);
        param_results = calc_param_ANOVA(mean_wind_data, [], (1:length(factor_levels > 1))+2, alpha);
        h_mean_amp(i) = param_results.h;
        
    end
    sim_time = toc;
    
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
        effect = 'none (null)';
        effect_description = 'all effects null';
    end
    
    %Print output
    fprintf('\n');
    fprintf('----------------------------------------------------------------------------------\n')
    fprintf('****SIMULATION SUMMARY****\n')
    fprintf('%s\n', time_stamp);
    fprintf('\nSimulation took %.2f minutes\n', sim_time/60);
    
    fprintf('\nUsing real EEG noise trials\n');
    fprintf('%d simulated experiments, %d permutations each, %d subjects\n', n_exp, n_perm, n_subs);
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
