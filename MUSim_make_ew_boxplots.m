%Make box plots for elementwise measures from MUSim simulations
%
%Author: Eric Fields
%Version Date: 17 June 2019

clearvars; close all;

%Make up some data
simulation_results.Fmax.ew_power = normrnd(0.2, 0.1, [1e4, 1]);
simulation_results.cluster01.ew_power = normrnd(0.4, 0.1, [1e4, 1]);
simulation_results.cluster05.ew_power = normrnd(0.6, 0.1, [1e4, 1]);
simulation_results.BH.ew_power = normrnd(0.5, 0.1, [1e4, 1]);
simulation_results.BY.ew_power = normrnd(0.1, 0.1, [1e4, 1]);
simulation_results.BKY.ew_power = normrnd(0.52, 0.1, [1e4, 1]);

methods = {'Fmax', 'cluster01', 'cluster05', 'BH', 'BY', 'BKY'};
names = {'Fmax', 'Cluster 0.01', 'Cluster 0.05', 'BH FDR', 'BY FDR', 'BKY FDR'};

%% Try 1

data = [];
for m = 1:length(methods)
    method = methods{m};
    data(:, m) = simulation_results.(method).ew_power;
    data(data(:, m)<0, m) = 0;
end

colors = [144 238 144; 0 0 128; 100 149 237; 255 0 0; 240 128 128; 178 34 34] ./ 255;

boxplot(data, string(names), ...
        'PlotStyle', 'traditional', ...
        'BoxStyle', 'filled', ...
        'Colors', colors);
%Make boxes wider
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',30); % Set width

%% Try 2

data = {};
for m = 1:length(methods)
    method = methods{m};
    data{m} = simulation_results.(method).ew_power;
end

py.matplotlib.pyplot.boxplot(data);