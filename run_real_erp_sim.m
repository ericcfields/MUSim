%Function to test Type I and Type II error rates for mass univariate approaches
%with data constructed from real EEG noise trials and real ERP effects
%
%Author: Eric Fields
%Version Date: 3 April 2019

function run_real_erp_sim(effect, time_wind, electrodes, factor_levels, mult_comp_method, n_exp, n_perm, save_results)

    
    %% ####################################################################
    %####################### PARAMETERS AND SET-UP ########################
    %######################################################################
    
    main_dir = MUSim_main_dir();
    
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; %#ok<ASGLU>
    close all;
    
    %% ~~~~~ SIMULATION PARAMETERS ~~~~~
    
    global VERBLEVEL
    VERBLEVEL = 0;
    
    %Load Emprob GND to use as basis for new GND
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
    
    %Save results?
    if n_exp > 1e3 && ~save_results
        resp = input('Are you sure you don''t want to save results?(y/n) ', 's');
        if strcmpi(resp, 'n')
            return;
        end
    elseif n_exp <= 1e3 && save_results
        resp = input('Are you sure you want to save the results?(y/n) ', 's');
        if strcmpi(resp, 'n')
            return;
        end
    end
    
    %Simulated data characteristics
    n_subs = 24;
    n_time_pts = length(GND.time_pts); %#ok<NODEF>
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
    n_param_time_pts = param_end_sample - param_start_sample + 1;
    
    %Other useful numbers
    n_electrodes = length(electrodes);
    n_conds      = length(bins);
    n_trials     = sum(cond_trials);
    

    %% ~~~~~ MAKE NEW GND ~~~~~
    
    %Update to match simulation parameters
    GND.grands = NaN(n_electrodes, n_time_pts, n_conds);
    GND.grands_stder = NaN(n_electrodes, n_time_pts, n_conds);
    GND.grands_t = NaN(n_electrodes, n_time_pts, n_conds);
    GND.sub_ct = ones(1, n_conds) * n_subs;
    GND.chanlocs = GND.chanlocs(electrodes);
    GND.bin_info = struct;
    for b = 1:length(bins)
        GND.bin_info(b).bindesc = char(64+b);
        GND.bin_info(b).condcode = 1;
    end
    GND.indiv_bin_ct = ones(n_subs, n_conds)*40;
    GND.indiv_bin_raw_ct = ones(n_subs, n_conds)*40;
    GND.indiv_erps = NaN(n_electrodes, n_time_pts, n_conds, n_subs);
    
 
    %% ####################################################################
    %#################### SIMULTATIONS ####################################
    %######################################################################
    
    %% ~~~~~ SET-UP ~~~~~
    
    %Find relevant time points for time window
    [~, start_sample] = min(abs(GND.time_pts - time_wind(1)));
    [~, end_sample  ] = min(abs(GND.time_pts - time_wind(2)));
    n_sample_time_pts = end_sample - start_sample + 1;
    
    %Pre-allocate variables
    F_int_nht            = NaN(n_exp, n_electrodes, n_sample_time_pts);
    sidak_int_nht        = NaN(n_exp, n_electrodes, n_sample_time_pts);
    param_uncorr_int_nht = NaN(n_exp, n_electrodes, n_sample_time_pts);
    param_ANOVA_nht      = NaN(1, n_exp);
    GND = repmat(GND, 1, n_exp);
    results = struct;
    results.clust05 = repmat(struct('h', NaN(n_electrodes, n_sample_time_pts), 'p', NaN(n_electrodes, n_sample_time_pts), ... 
                             'F_obs', NaN(n_electrodes, n_sample_time_pts),'df', NaN(1, length(factor_levels)), ... 
                              'clust_info', struct('null_test', [], 'pval', [], 'clust_mass', [], 'clust_ids', []), ...
                              'estimated_alpha', NaN, 'exact_test', NaN), ...
                              1, n_exp);
    results.clust01 = results.clust05;
    results.Fmax = repmat(struct('h', NaN(n_electrodes, n_sample_time_pts), 'p', NaN(n_electrodes, n_sample_time_pts), ... 
                          'F_obs', NaN(n_electrodes, n_sample_time_pts),'Fmax_crit', NaN, 'df', NaN(1, length(factor_levels)), ...
                          'estimated_alpha', NaN, 'exact_test', NaN), ...
                          1, n_exp);
    results.bh = repmat(struct('h', NaN(n_electrodes, n_time_pts), 'p', NaN(n_electrodes, n_time_pts), ... 
                        'F_obs', NaN(n_electrodes, n_time_pts), 'F_crit', NaN, 'df', NaN(1, 2)), ...
                        1, n_exp);
    results.by  = results.bh;
    results.bky = results.bh;

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
    for i = 1:n_exp
        
        %% ~~~~~ GENERATE SIMULATED DATA ~~~~~
        
        %Randomly select trials from each subject
        for s = 1:n_subs
            GND(i).indiv_erps(:,:,:,s) = mean(reshape(all_noise_trials{sub_sample(i,s),2}(electrodes, :, randsample(1:size(all_noise_trials{sub_sample(i,s),2},3), n_trials)), n_electrodes, n_time_pts, n_conds, []), 4); %#ok<PFIIN,PFOUS>
        end
        
        %Add effects if any
        if effect
            %Generate subject effects for this experiment
            sub_effects = normrnd(sub_param(1), sub_param(2), n_subs, n_conds);
            for s = 1:n_subs
                for b = 1:n_conds
                    GND(i).indiv_erps(:, :, b, s) = GND(i).indiv_erps(:, :, b, s)*error_mult(b) + effects_data(:, :, b)*sub_effects(s, b);
                end
            end
        end

%         %Use data to populate other GND fields
%         GND.grands = mean(GND.indiv_erps, 4);
%         GND.grands_stder = std(GND.indiv_erps, 0, 4) ./ sqrt(size(GND.indiv_erps, 4));
%         GND.grands_t = GND.grands ./ GND.grands_stder;
        
        %% ~~~~~ RUN STATS ~~~~~
        
        %ANOVA
        data = GND(i).indiv_erps(:, start_sample:end_sample, bins, :);
        data = reshape(data,[n_electrodes, n_sample_time_pts, factor_levels, n_subs]);
        results(i).Fmax    = calc_Fmax(data, [], (1:length(factor_levels > 1))+2, n_perm, 0.05); %#ok<PFOUS>
        results(i).clust05 = calc_Fclust(data, [], (1:length(factor_levels > 1))+2, n_perm, 0.05, chan_hood, thresh_p);
        results(i).clust01 = calc_Fclust(data, [], (1:length(factor_levels > 1))+2, n_perm, 0.01, chan_hood, thresh_p);
        results(i).bh      = calc_param_ANOVA(data, [], (1:length(factor_levels > 1))+2, 0.05, 'bh');
        results(i).by      = calc_param_ANOVA(data, [], (1:length(factor_levels > 1))+2, 0.05, 'by');
        results(i).bky     = calc_param_ANOVA(data, [], (1:length(factor_levels > 1))+2, 0.05, 'bky');
        
        %Parametric ANOVA
        if isequal(electrodes, 1:32)
            mean_wind_data = mean(mean(GND(i).indiv_erps(param_electrodes, param_start_sample:param_end_sample, bins, :), 1), 2);
            mean_wind_data = reshape(mean_wind_data, [size(mean_wind_data, 1), size(mean_wind_data, 2), factor_levels, size(mean_wind_data, 4)]);
        elseif isequal(electrodes, param_electrodes)
            mean_wind_data = mean(mean(GND(i).indiv_erps(:, param_start_sample:param_end_sample, bins, :), 1), 2);
            mean_wind_data = reshape(mean_wind_data, [size(mean_wind_data, 1), size(mean_wind_data, 2), factor_levels, size(mean_wind_data, 4)]);
        else
            error('Can''t determine electrodes for parametric test');
        end
        param_results = calc_param_ANOVA(mean_wind_data, [], (1:length(factor_levels > 1))+2, 0.05);
        param_ANOVA_nht(i) = param_results.h;
                             
        %Store hypothesis test results
        F_int_nht(i, :, :)            = results(i).(mult_comp_method).h;
        sidak_int_nht(i, :, :)        = results(i).(mult_comp_method).F_obs > finv(.95^(1/(n_electrodes * n_sample_time_pts)), prod(factor_levels-1), prod(factor_levels-1) * (n_subs - 1));
        param_uncorr_int_nht(i, :, :) = results(i).(mult_comp_method).F_obs > finv(.95,                                        prod(factor_levels-1), prod(factor_levels-1) * (n_subs - 1));

    end
    toc
    
    
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
    
    %Temporal extent
    sig_exp = F_int_nht(any(F_int_nht(:, :)'), :, :);
    temporal_extent = sum((any(sig_exp, 2)), 3);

    
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
    if strcmpi(mult_comp_method, 'clust')
        fprintf('\nThreshold p = %.3f', thresh_p);
    end
    fprintf('\nMultiple comparison method: %s\n\n', mult_comp_method);
    
    fprintf('CORRECTED F-TEST RESULTS\n')
    fprintf('Rejection rate all time points (individually) = %f\n', mean(F_int_nht(:)));
    fprintf('Family-wise rejection rate = %f\n\n', mean(any(reshape(F_int_nht, n_exp, n_electrodes*n_sample_time_pts)')));
    
    fprintf('MEAN WINDOW/REGION PARAMETRIC F-TEST RESULTS\n')
    fprintf('Electrodes: ');
    fprintf([sprintf('%d, ', param_electrodes(1:end-1)), num2str(param_electrodes(end))]);
    fprintf('\nTime window: %d - %d\n', param_time_wind(1), param_time_wind(2));
    fprintf('Rejection rate = %f\n\n', mean(param_ANOVA_nht));
    
    fprintf('DUNN-SIDAK CORRECTED PARAMETRIC F-TEST RESULTS\n');
    fprintf('Rejection rate all time points (individually) = %f\n', mean(sidak_int_nht(:)));
    fprintf('Family-wise rejection rate = %f\n\n', mean(any(reshape(sidak_int_nht, n_exp, n_electrodes*n_sample_time_pts)')));
    
    fprintf('UNCORRECTED PARAMETRIC F-TEST RESULTS\n');
    fprintf('Rejection rate all time points (individually) = %f\n', mean(param_uncorr_int_nht(:)));
    fprintf('Family-wise rejection error rate = %f\n\n', mean(any(reshape(param_uncorr_int_nht, n_exp, n_electrodes*n_sample_time_pts)')));
    
    fprintf('TEMPORAL EXTENT\n');
    fprintf('Mean length of effects was %.0f ms\n', mean(temporal_extent)*(1000/128));
    fprintf('Median length of effects was %.0f ms\n\n', median(temporal_extent)*(1000/128));
    
    fprintf('----------------------------------------------------------------------------------\n\n')
    
    diary off
    
end
