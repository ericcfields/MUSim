%Read simulation results in csv files
%
%Author: Eric Fields
%Version Date: 20 June 2019

main_dir = MUSim_main_dir();

results_mats = get_files(fullfile(main_dir, 'results'), 'simulation_results.mat');

components = {'namesP300_reduced', 'P300'; ...
              'NonCon_N400_reduced', 'N400'; ...
              'simulated_focal_effect', 'P1'};

for i = 1
    
    mat_file = results_mats{i};
    
    %Get effect and time window
    effect_name = components{cellfun(@(x) contains(mat_file, x), components(:,1)), 2};
    [s_idx, e_idx] = regexp(mat_file, '_\d*-\d*_');
    time_wind = mat_file((s_idx+1):(e_idx-1));
    
    %Load simulation results
    load(fullfile(main_dir, 'results', mat_file));
    methods = fieldnames(simulation_results);
    
    %Elementwise Power
    ew_power = table('Size', [1e4, length(methods)-1], ... 
                     'VariableTypes', repmat({'doublenan'}, [1, length(methods)-1]), ...
                     'VariableNames', methods(2:end));
    for m = 2:length(methods)
        method = methods{m};
        ew_power{1:length(simulation_results.(method).ew_power), methods{m}} = simulation_results.(method).ew_power;
    end
    
end
