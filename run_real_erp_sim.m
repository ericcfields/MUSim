%Test Type I and Type II error rates for mass univariate approaches
%with data constructed from real EEG noise trials and real ERP effects
%
%Author: Eric Fields
%Version Date: 9 April 2019

function run_real_erp_sim(noise, effect, time_wind, electrodes, factor_levels, n_exp, n_perm, n_subs, cond_trials, error_mult, ind_var_factor, alpha, output_file)

    
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
        [F_obs, F_dist, df_effect, df_res] = perm_rbANOVA(data, 3, n_perm);
        
        %Fmax
        h_Fmax(i, :, :) = Fmax_corr(F_obs, F_dist, alpha);
        thresh_F = finv(1-0.05, df_effect, df_res);
        
        %cluster
        h_clust05(i, :, :) = Fclust_corr(F_obs, F_dist, alpha, chan_hood, thresh_F);
        thresh_F = finv(1-0.01, df_effect, df_res);
        h_clust01(i, :, :) = Fclust_corr(F_obs, F_dist, alpha, chan_hood, thresh_F);
        
        %FDR
        p_uncorrected(i, :, :) = 1 - fcdf(F_obs, df_effect, df_res);
        h_uncorrected(i, :, :) = p_uncorrected(i, :, :) <= alpha;
        h_bh(i, :, :)  = fdr_bh(p_uncorrected(i, :, :), alpha, 'pdep', 'no');
        h_by(i, :, :)  = fdr_bh(p_uncorrected(i, :, :), alpha, 'dep', 'no');
        h_bky(i, :, :) = fdr_bky(p_uncorrected(i, :, :), alpha, 'no');
        
        %Mean amplitude ANOVA
        mean_wind_data = mean(mean(data, 1), 2);
        mean_wind_data = reshape(mean_wind_data, [size(mean_wind_data, 1), size(mean_wind_data, 2), factor_levels, size(mean_wind_data, 4)]);
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
    
    %Save results
    time_stamp = sprintf('%s_%s', datestr(datetime('now'),'ddmmmyy'), datestr(datetime('now'),'HHMM'));
    if output_file
        diary(output_file);
    end
    
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
    fprintf('Rejection rate =\t%.3f\n', mean(h_mean_amp));
    
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
    summarize_results(effect_loc, h_uncorrected);
    
    fprintf('\nSIDAK RESULTS\n');
    summarize_results(effect_loc, h_sidak);
    
    fprintf('\nFMAX RESULTS\n');
    summarize_results(effect_loc, h_Fmax);
    
    fprintf('\nCLUSTER 0.05 RESULTS\n');
    summarize_results(effect_loc, h_clust05);
    
    fprintf('\nCLUSTER 0.01 RESULTS\n');
    summarize_results(effect_loc, h_clust01);
    
    fprintf('\nBH FDR RESULTS\n');
    summarize_results(effect_loc, h_bh);
    
    fprintf('\nBY FDR RESULTS\n');
    summarize_results(effect_loc, h_by);
    
    fprintf('\nBKY FDR RESULTS\n');
    summarize_results(effect_loc, h_bky);
    
    fprintf('\n----------------------------------------------------------------------------------\n\n')
    
    diary off
    
end

function summarize_results(effect_loc, nht)

    [n_perm, ~, n_time_pts] = size(nht);

    %Get simulated experiments that found a significant result
    sig_studies = any(any(nht, 2), 3);
    
    %Collapse across electrodes
    nht_t = reshape(any(nht, 2), [n_perm, n_time_pts]);
    
    %Get null hypothesis test at locations with and without real effect
    %separately
    nht_effect = nht_t(:, effect_loc);
    nht_null   = nht_t(:, ~effect_loc);
    
    %Report family-wiise rejection rate
    fprintf('-- Family-wise rejection rates --\n');
    fprintf('Family-wise rejection rate across time points with effect (familywise power) =\t%.3f\n',                 mean(any(nht_effect, 2)));
    fprintf('Family-wise rejection rate across time points with null effect (familywise Type I error) =\t%.3f\n',     mean(any(nht_null, 2)));
    
    %Report element-wise rejection rate within studies with significant
    %results
    fprintf('-- Element-wise rejection rates --\n');
    fprintf('Mean rejection rate at individual time points with effect (element-wise power) =\t%.3f\n',               mean(mean(nht_effect(sig_studies, :))));
    fprintf('Median rejection rate at individual time points with effect (element-wise power) =\t%.3f\n',             median(mean(nht_effect(sig_studies, :), 2)));
    fprintf('Mean rejection rate at individual time points with null effect (element-wise Type I error) =\t%.3f\n',   mean(mean(nht_null(sig_studies, :))));
    fprintf('Median rejection rate at individual time points with null effect (element-wise Type I error) =\t%.3f\n', median(mean(nht_null(sig_studies, :), 2)));
    fprintf('Mean element-wise false discovery rate =\t%.3f\n',                                                       mean(sum(nht_null, 2) ./ (sum(nht_null, 2) + sum(nht_effect, 2))));
    fprintf('Median element-wise false discovery rate =\t%.3f\n',                                                     median(sum(nht_null, 2) ./ (sum(nht_null, 2) + sum(nht_effect, 2))));
    
end
