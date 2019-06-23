%Read simulation results in csv files
%
%Author: Eric Fields
%Version Date: 21 June 2019

main_dir = MUSim_main_dir();

%% Familywise measures
%Run Python script creating Familywise measure csvs

%Add location of Python script to Python path
py_addpath(main_dir);

%Run script
py.MUSim_make_csv.main();


%% Elementwise measurse

results_mats = get_files(fullfile(main_dir, 'results'), 'simulation_results.mat');

components = {'namesP300_reduced', 'P300'; ...
              'NonCon_N400_reduced', 'N400'; ...
              'simulated_focal_effect', 'P1'};

for i = 1:length(results_mats)
    
    mat_file = results_mats{i};
    
    %Get effect and time window
    effect_name = components{cellfun(@(x) contains(mat_file, x), components(:,1)), 2};
    [s_idx, e_idx] = regexp(mat_file, '_\d*-\d*_');
    time_wind = mat_file((s_idx+1):(e_idx-1));
    
    %Load simulation results
    load(fullfile(main_dir, 'results', mat_file));
    methods = fieldnames(simulation_results);
    methods = methods(9:end);
    
    %Elementwise Power
    ew_power = table('Size', [1e4, length(methods)-1], ... 
                     'VariableTypes', repmat({'doublenan'}, [1, length(methods)-1]), ...
                     'VariableNames', methods(2:end));
    for m = 2:length(methods)
        method = methods{m};
        ew_power{:, methods{m}} = simulation_results.(method).ew_power;
        output_file = fullfile(main_dir, 'results', sprintf('MUSim_Power_EW_power_%s_%s.csv', effect_name, time_wind));
        writetable(ew_power, output_file);
    end
    
    %Elementwise FDR
    ew_FDR = table('Size', [1e4, length(methods)-1], ... 
                     'VariableTypes', repmat({'doublenan'}, [1, length(methods)-1]), ...
                     'VariableNames', methods(2:end));
    for m = 2:length(methods)
        method = methods{m};
        ew_FDR{:, methods{m}} = simulation_results.(method).ew_FDR;
        output_file = fullfile(main_dir, 'results', sprintf('MUSim_Power_EW_FDR_%s_%s.csv', effect_name, time_wind));
        writetable(ew_FDR, output_file);
    end
    
    %Elementwise Type I
    ew_TypeI = table('Size', [1e4, length(methods)-1], ... 
                     'VariableTypes', repmat({'doublenan'}, [1, length(methods)-1]), ...
                     'VariableNames', methods(2:end));
    for m = 2:length(methods)
        method = methods{m};
        ew_TypeI{:, methods{m}} = simulation_results.(method).ew_TypeI;
        output_file = fullfile(main_dir, 'results', sprintf('MUSim_Power_EW_TypeI_%s_%s.csv', effect_name, time_wind));
        writetable(ew_TypeI, output_file);
    end
    
    %Onset
    onset = table('Size', [1e4, length(methods)-1], ... 
                     'VariableTypes', repmat({'doublenan'}, [1, length(methods)-1]), ...
                     'VariableNames', methods(2:end));
    for m = 2:length(methods)
        method = methods{m};
        onset{:, methods{m}} = simulation_results.(method).onset;
        output_file = fullfile(main_dir, 'results', sprintf('MUSim_Power_EW_onset_%s_%s.csv', effect_name, time_wind));
        writetable(onset, output_file);
    end
    
    %Offset
    offset = table('Size', [1e4, length(methods)-1], ... 
                     'VariableTypes', repmat({'doublenan'}, [1, length(methods)-1]), ...
                     'VariableNames', methods(2:end));
    for m = 2:length(methods)
        method = methods{m};
        offset{:, methods{m}} = simulation_results.(method).offset;
        output_file = fullfile(main_dir, 'results', sprintf('MUSim_Power_EW_offset_%s_%s.csv', effect_name, time_wind));
        writetable(offset, output_file);
    end
    
end
