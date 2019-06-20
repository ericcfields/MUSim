%Test Type I and Type II error rates for mass univariate approaches
%with data constructed from real EEG noise trials and real ERP effects
%
%Author: Eric Fields
%Version Date: 20 June 2019
%
%Copyright (c) 2019, Eric C. Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.
%This software is provided "as is" and any express or implied warranties are disclaimed. 

function run_real_erp_sim(noise, effect, time_wind, electrodes, factor_levels, dims, n_exp, n_perm, n_subs, cond_trials, error_mult, ind_var_factor, alpha, text_output, mat_output)

    
    %% ####################################################################
    %####################### PARAMETERS AND SET-UP ########################
    %######################################################################
    
    %Suppress command line output of stats functions
    global VERBLEVEL
    VERBLEVEL = 0;
    
    %Add EEGLAB, MUT, and FMUT functions to path
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; %#ok<ASGLU>
    close all;
    
    %Load EEG noise trials
    load(noise, 'noise_trials');
    noise_trials = noise_trials; %#ok<ASGSL>
    
    %Load effect data
    if ~strcmpi(effect, 'null')
        load(effect, 'effects_data', 'effect_description');
        if ~isequal(factor_levels, size(effects_data, 3))
            error('factor_levels input doesn''t match effects data');
        end
        effects_data = effects_data(electrodes, :, :);
    else
        effects_data = [];
    end

    %Some key numbers
    n_conds      = prod(factor_levels);
    n_time_pts   = length(noise_trials.times);
    n_electrodes = length(electrodes);
    
    %Multiplier applied to standard deviation of each bin
    if isscalar(error_mult)
        error_mult  = ones(n_conds, 1) * error_mult; 
    elseif length(error_mult) ~= prod(factor_levels)
        error('error_mult input must be a scalar or length must match the number of conditions');
    end
    
    %Number of trials in each condition
    if isscalar(cond_trials)
        cond_trials = ones(n_conds, 1) * cond_trials; %Number of trials in each bin
    elseif length(cond_trials) ~= prod(factor_levels)
        error('n_trials input must be a scalar or length must match the number of conditions');
    end
    %total number of trials
    n_trials = sum(cond_trials);
    
    %Channel neighbor information for cluster tests
    chan_hood = 75;
    chan_hood = spatial_neighbors(noise_trials.chanlocs(electrodes), chan_hood, []);
    
    %Find index of start and stop time point
    [~, start_sample] = min(abs(noise_trials.times - time_wind(1)));
    [~, end_sample  ] = min(abs(noise_trials.times - time_wind(2)));
    n_sample_time_pts = end_sample - start_sample + 1;
    
    %Pre-allocate output variables
    h_uncorrected  = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_sidak        = NaN(n_exp, n_electrodes, n_sample_time_pts); %#ok<NASGU>
    h_Fmax         = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_clust05      = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_clust01      = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_bh           = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_by           = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_bky          = NaN(n_exp, n_electrodes, n_sample_time_pts);
    h_mean_amp = NaN(n_exp, 1);
    p_uncorrected = NaN(n_exp, n_electrodes, n_sample_time_pts);
    
    
    %% ####################################################################
    %#################### SIMULTATIONS ####################################
    %######################################################################
    
    tic
    %Conduct n_exp simulated experiments
    parfor i = 1:n_exp
        
        %~~~~~~~~~~~~~~~~~~~~~~~ GENERATE SIMULATED DATA ~~~~~~~~~~~~~~~~~~~~~~~
        
        %Select random subset of subjects
        sub_sample = randsample(1:length(noise_trials.data), n_subs); %#ok<PFBNS>
        
        %Randomly select trials from each subject
        sim_data = NaN(n_electrodes, n_time_pts, n_conds, n_subs);
        for s = 1:n_subs
            sub_data = noise_trials.data{sub_sample(s)};
            sim_data(:,:,:,s) = mean(reshape(sub_data(electrodes, :, randsample(1:size(sub_data,3), n_trials)), n_electrodes, n_time_pts, n_conds, []), 4);
        end
        
        %Add effects if any
        if ~strcmpi(effect, 'null')
            %Generate subject effects for this experiment
            sub_effects = normrnd(1, ind_var_factor, n_subs, n_conds);
            for s = 1:n_subs
                for b = 1:n_conds
                    sim_data(:, :, b, s) = sim_data(:, :, b, s)*error_mult(b) + effects_data(:, :, b)*sub_effects(s, b); %#ok<PFBNS>
                end
            end
        end
        
        % ~~~~~~~~~~~~~~~~~~~~~~~ RUN STATS ~~~~~~~~~~~~~~~~~~~~~~~
        
        %Get analysis data
        data = sim_data(:, start_sample:end_sample, :, :);
        data = reshape(data,[n_electrodes, n_sample_time_pts, factor_levels, n_subs]);
        
        %Run ANOVA
        [F_obs, F_dist, df_effect, df_res] = perm_rbANOVA(data, dims, n_perm);
        
        %Uncorrected
        p_uncorrected(i, :, :) = 1 - fcdf(F_obs, df_effect, df_res);
        h_uncorrected(i, :, :) = p_uncorrected(i, :, :) <= alpha;
        
        %Fmax
        h_Fmax(i, :, :) = Fmax_corr(F_obs, F_dist, alpha);
        
        %cluster
        thresh_F = finv(1-0.05, df_effect, df_res);
        h_clust05(i, :, :) = Fclust_corr(F_obs, F_dist, alpha, chan_hood, thresh_F);
        thresh_F = finv(1-0.01, df_effect, df_res);
        h_clust01(i, :, :) = Fclust_corr(F_obs, F_dist, alpha, chan_hood, thresh_F);
        
        %FDR
        h_bh(i, :, :)  = fdr_bh(p_uncorrected(i, :, :), alpha, 'pdep', 'no');
        h_by(i, :, :)  = fdr_bh(p_uncorrected(i, :, :), alpha, 'dep', 'no');
        h_bky(i, :, :) = fdr_bky(p_uncorrected(i, :, :), alpha, 'no');
        
        %Mean amplitude ANOVA
        mean_wind_data = mean(mean(data, 1), 2);
        mean_wind_data = reshape(mean_wind_data, [size(mean_wind_data, 1), size(mean_wind_data, 2), factor_levels, size(mean_wind_data, ndims(mean_wind_data))]);
        param_results = calc_param_ANOVA(mean_wind_data, [], (1:length(factor_levels > 1))+2, alpha);
        h_mean_amp(i) = param_results.h;
        
    end
    sim_time = toc;
    
    %Calcluate Dunn-Sidak corrected results
    p_sidak = 1-(1-p_uncorrected).^(n_electrodes*n_time_pts);
    p_sidak(p_sidak>1) = 1;
    h_sidak = p_sidak <= alpha;
    
    
    %% ####################################################################
    %################# SAVE AND REPORT RESULTS ############################
    %######################################################################
    
    %Save results to text file
    time_stamp = sprintf('%s_%s', datestr(datetime('now'),'ddmmmyy'), datestr(datetime('now'),'HHMM'));
    if text_output
        diary(text_output);
    end
    
    %Save results to .mat file
    simulation_results = struct;
    simulation_results.effect = effect;
    simulation_results.effect_description = effect_description;
    simulation_results.time_window = time_wind;
    simulation_results.electrodes = electrodes;
    simulation_results.n_experiments = n_exp;
    simulation_results.n_permutations = n_perm;
    simulation_results.n_subjects = n_subs;
    simulation_results.n_trials = cond_trials;
    
    if strcmpi(effect, 'null')
        effect_description = '';
    end
    
    %Print output
    fprintf('\n');
    fprintf('----------------------------------------------------------------------------------\n')
    fprintf('****SIMULATION SUMMARY****\n')
    fprintf('%s\n', time_stamp);
    fprintf('\nSimulation took\t%.2f minutes\n', sim_time/60);
    
    fprintf('\nUsing real EEG noise trials\n');
    fprintf('Simulated experiments =\t%d\n', n_exp);
    fprintf('Permutations =\t%d\n', n_perm);
    fprintf('Sample size =\t%d\n', n_subs);
    fprintf('\nEffect:\t%s', effect);
    fprintf('\nEffect description:\t%s\n', effect_description);
    fprintf('Time window:\t%d - %d\n', time_wind(1), time_wind(2));
    fprintf('Electrodes:\t');
    fprintf([sprintf('%d, ', electrodes(1:end-1)), num2str(electrodes(end))]);
    fprintf('\nFactor levels:\t');
    fprintf('%d  ', factor_levels);
    fprintf('\nError multiplier =\t');
    fprintf('%.1f  ', error_mult);
    fprintf('\nTrials =\t');
    fprintf('%d  ', cond_trials);
    fprintf('\n');
    
    fprintf('\nMEAN WINDOW/REGION PARAMETRIC F-TEST RESULTS\n')
    simulation_results.mean_amp.rej_rate = mean(h_mean_amp);
    fprintf('Rejection rate =\t%.3f\n', simulation_results.mean_amp.rej_rate);  
    
    %Find time points with effect
    effect_loc = false(1, end_sample-start_sample+1);
    if ~strcmpi(effect, 'null')
        i = 0;
        for t = start_sample:end_sample
            i = i + 1;
            effect_loc(i) = ~all(all(effects_data(:, t, :) == effects_data(:, t, 1)));
        end
    end
    
    fprintf('\nUNCORRECTED RESULTS\n');
    [fw_power, fw_TypeI, fw_total_miss, fw_FDR, ew_power, ew_TypeI, ew_FDR, onset_time, offset_time] = summarize_results(effect_loc, h_uncorrected, noise_trials.times(start_sample:end_sample));
    simulation_results.uncorrected.fw_power = fw_power;
    simulation_results.uncorrected.fw_TypeI = fw_TypeI;
    simulation_results.uncorrected.fw_total_miss = fw_total_miss;
    simulation_results.uncorrected.fw_FDR = fw_FDR;
    simulation_results.uncorrected.ew_power = ew_power;
    simulation_results.uncorrected.ew_TypeI = ew_TypeI;
    simulation_results.uncorrected.ew_FDR = ew_FDR;
    simulation_results.uncorrected.onset = onset_time;
    simulation_results.uncorrected.offset = offset_time;
    
    fprintf('\nSIDAK RESULTS\n');
    [fw_power, fw_TypeI, fw_total_miss, fw_FDR, ew_power, ew_TypeI, ew_FDR, onset_time, offset_time] = summarize_results(effect_loc, h_sidak, noise_trials.times(start_sample:end_sample));
    simulation_results.sidak.fw_power = fw_power;
    simulation_results.sidak.fw_TypeI = fw_TypeI;
    simulation_results.sidak.fw_total_miss = fw_total_miss;
    simulation_results.sidak.fw_FDR = fw_FDR;
    simulation_results.sidak.ew_power = ew_power;
    simulation_results.sidak.ew_TypeI = ew_TypeI;
    simulation_results.sidak.ew_FDR = ew_FDR;
    simulation_results.sidak.onset = onset_time;
    simulation_results.sidak.offset = offset_time;
    
    fprintf('\nFMAX RESULTS\n');
    [fw_power, fw_TypeI, fw_total_miss, fw_FDR, ew_power, ew_TypeI, ew_FDR, onset_time, offset_time] = summarize_results(effect_loc, h_Fmax, noise_trials.times(start_sample:end_sample));
    simulation_results.Fmax.fw_power = fw_power;
    simulation_results.Fmax.fw_TypeI = fw_TypeI;
    simulation_results.Fmax.fw_total_miss = fw_total_miss;
    simulation_results.Fmax.fw_FDR = fw_FDR;
    simulation_results.Fmax.ew_power = ew_power;
    simulation_results.Fmax.ew_TypeI = ew_TypeI;
    simulation_results.Fmax.ew_FDR = ew_FDR;
    simulation_results.Fmax.onset = onset_time;
    simulation_results.Fmax.offset = offset_time;
    
    fprintf('\nCLUSTER 0.05 RESULTS\n');
    [fw_power, fw_TypeI, fw_total_miss, fw_FDR, ew_power, ew_TypeI, ew_FDR, onset_time, offset_time] = summarize_results(effect_loc, h_clust05, noise_trials.times(start_sample:end_sample));
    simulation_results.cluster05.fw_power = fw_power;
    simulation_results.cluster05.fw_TypeI = fw_TypeI;
    simulation_results.cluster05.fw_total_miss = fw_total_miss;
    simulation_results.cluster05.fw_FDR = fw_FDR;
    simulation_results.cluster05.ew_power = ew_power;
    simulation_results.cluster05.ew_TypeI = ew_TypeI;
    simulation_results.cluster05.ew_FDR = ew_FDR;
    simulation_results.cluster05.onset = onset_time;
    simulation_results.cluster05.offset = offset_time;
    
    fprintf('\nCLUSTER 0.01 RESULTS\n');
    [fw_power, fw_TypeI, fw_total_miss, fw_FDR, ew_power, ew_TypeI, ew_FDR, onset_time, offset_time] = summarize_results(effect_loc, h_clust01, noise_trials.times(start_sample:end_sample));
    simulation_results.cluster01.fw_power = fw_power;
    simulation_results.cluster01.fw_TypeI = fw_TypeI;
    simulation_results.cluster01.fw_total_miss = fw_total_miss;
    simulation_results.cluster01.fw_FDR = fw_FDR;
    simulation_results.cluster01.ew_power = ew_power;
    simulation_results.cluster01.ew_TypeI = ew_TypeI;
    simulation_results.cluster01.ew_FDR = ew_FDR;
    simulation_results.cluster01.onset = onset_time;
    simulation_results.cluster01.offset = offset_time;
    
    fprintf('\nBH FDR RESULTS\n');
    [fw_power, fw_TypeI, fw_total_miss, fw_FDR, ew_power, ew_TypeI, ew_FDR, onset_time, offset_time] = summarize_results(effect_loc, h_bh, noise_trials.times(start_sample:end_sample));
    simulation_results.BH.fw_power = fw_power;
    simulation_results.BH.fw_TypeI = fw_TypeI;
    simulation_results.BH.fw_total_miss = fw_total_miss;
    simulation_results.BH.fw_FDR = fw_FDR;
    simulation_results.BH.ew_power = ew_power;
    simulation_results.BH.ew_TypeI = ew_TypeI;
    simulation_results.BH.ew_FDR = ew_FDR;
    simulation_results.BH.onset = onset_time;
    simulation_results.BH.offset = offset_time;
    
    fprintf('\nBY FDR RESULTS\n');
    [fw_power, fw_TypeI, fw_total_miss, fw_FDR, ew_power, ew_TypeI, ew_FDR, onset_time, offset_time] = summarize_results(effect_loc, h_by, noise_trials.times(start_sample:end_sample));
    simulation_results.BY.fw_power = fw_power;
    simulation_results.BY.fw_TypeI = fw_TypeI;
    simulation_results.BY.fw_total_miss = fw_total_miss;
    simulation_results.BY.fw_FDR = fw_FDR;
    simulation_results.BY.ew_power = ew_power;
    simulation_results.BY.ew_TypeI = ew_TypeI;
    simulation_results.BY.ew_FDR = ew_FDR;
    simulation_results.BY.onset = onset_time;
    simulation_results.BY.offset = offset_time;
    
    
    fprintf('\nBKY FDR RESULTS\n');
    [fw_power, fw_TypeI, fw_total_miss, fw_FDR, ew_power, ew_TypeI, ew_FDR, onset_time, offset_time] = summarize_results(effect_loc, h_bky, noise_trials.times(start_sample:end_sample));
    simulation_results.BKY.fw_power = fw_power;
    simulation_results.BKY.fw_TypeI = fw_TypeI;
    simulation_results.BKY.fw_total_miss = fw_total_miss;
    simulation_results.BKY.fw_FDR = fw_FDR;
    simulation_results.BKY.ew_power = ew_power;
    simulation_results.BKY.ew_TypeI = ew_TypeI;
    simulation_results.BKY.ew_FDR = ew_FDR;
    simulation_results.BKY.onset = onset_time;
    simulation_results.BKY.offset = offset_time;
    
    fprintf('\n----------------------------------------------------------------------------------\n\n')
    
    diary off
    
    %Save results struct
    save(mat_output, 'simulation_results');
    
end

function [fw_power, fw_TypeI, fw_total_miss, fw_FDR, ew_power, ew_TypeI, ew_FDR, onset_time, offset_time] = summarize_results(effect_loc, nht, time_ids)

    [n_exp, ~, n_time_pts] = size(nht);

    %Get simulated experiments that found a significant result
    sig_studies = any(any(nht, 2), 3);
    
    %Collapse across electrodes
    nht_t = reshape(any(nht, 2), [n_exp, n_time_pts]);
    
    %Get null hypothesis test at locations with and without real effect
    %separately
    nht_effect = nht_t(:, effect_loc);
    nht_null   = nht_t(:, ~effect_loc);
    
    %Get the same info for the subset of studies that found significant
    %results
    nht_effect_sig = nht_effect(sig_studies, :);
    nht_null_sig   = nht_null(sig_studies, :);
    
    %Report family-wiise rejection rates
    fprintf('-- Family-wise rejection rates --\n');
    nht_effect_fw = any(nht_effect, 2);
    nht_null_fw   = any(nht_null, 2);
    fw_power      = mean(nht_effect_fw);
    fw_TypeI      = mean(nht_null_fw);
    fw_total_miss = mean(~nht_effect_fw & nht_null_fw);
    fw_FDR        = sum(nht_null_fw)/sum(nht_effect_fw | nht_null_fw);
    fprintf('Family-wise rejection rate across time points with effect (familywise power) =\t%.3f\n', fw_power);
    fprintf('Family-wise rejection rate across time points with null effect (familywise Type I error) =\t%.3f\n', fw_TypeI);
    fprintf('Total miss rate (only null time points rejected) =\t%.3f\n', fw_total_miss);
    fprintf('Familywise FDR (proportion of sig studies that include false positive time point) =\t%.3f\n', fw_FDR);
    
    %Report element-wise rejection rates within studies with significant
    %results
    fprintf('-- Element-wise rejection rates --\n');
    ew_power = NaN(n_exp, 1);
    ew_power(sig_studies) = mean(nht_effect_sig, 2);
    ew_TypeI = NaN(n_exp, 1);
    ew_TypeI(sig_studies) = mean(nht_null_sig, 2);
    ew_FDR = NaN(n_exp, 1);
    ew_FDR(sig_studies)   = sum(nht_null_sig, 2) ./ (sum(nht_null_sig, 2) + sum(nht_effect_sig, 2));
    fprintf('Mean rejection rate at individual time points with effect (element-wise power) =\t%.3f\n',               nanmean(ew_power));
    fprintf('Median rejection rate at individual time points with effect (element-wise power) =\t%.3f\n',             nanmedian(ew_power));
    fprintf('Mean rejection rate at individual time points with null effect (element-wise Type I error) =\t%.3f\n',   nanmean(ew_TypeI));
    fprintf('Median rejection rate at individual time points with null effect (element-wise Type I error) =\t%.3f\n', nanmedian(ew_TypeI));
    fprintf('Mean element-wise false discovery rate =\t%.3f\n',                                                       nanmean(ew_FDR));
    fprintf('Median element-wise false discovery rate =\t%.3f\n',                                                     nanmedian(ew_FDR));
    
    %Get earliest and latest time points
    fprintf('-- Onset and Offset Times --\n');
    if any(sig_studies)
        
        sig_onset_time = NaN(sum(sig_studies), 1);
        sig_offset_time = NaN(sum(sig_studies), 1);
        j = 0;
        for i = 1:n_exp
            if sig_studies(i)
                j = j + 1;
                sig_onset_time(j) = find(nht_t(i, :), 1);
                sig_offset_time(j) = find(nht_t(i, :), 1, 'last');
            end
        end
        
        onset_time = NaN(n_exp, 1);
        onset_time(sig_studies) = sig_onset_time;
        offset_time = NaN(n_exp, 1);
        offset_time(sig_studies) = sig_offset_time;
        
        fprintf('Mean onset time =\t%d\n', time_ids(round(nanmean(onset_time))));
        fprintf('Median onset time =\t%d\n', time_ids(round(nanmedian(onset_time))));
        fprintf('Mean offset time =\t%d\n', time_ids(round(nanmean(offset_time))));
        fprintf('Median offset time =\t%d\n', time_ids(round(nanmedian(offset_time))));
        
    else
        
        fprintf('Mean onset time =\t%d\n', NaN);
        fprintf('Median onset time =\t%d\n', NaN);
        fprintf('Mean offset time =\t%d\n', NaN);
        fprintf('Median offset time =\t%d\n', NaN);
        
    end
    
end
