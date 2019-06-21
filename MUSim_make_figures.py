# -*- coding: utf-8 -*-
"""
Make figures for MUSim paper

AUTHOR: Eric Fields
VERSION DATE: 15 June 2019
"""

import os
from os.path import join
import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportion_confint
import matplotlib.pyplot as plt

results_dir = r'C:\Users\ecfne\Documents\Eric\Research\Stats Simulations\MUSim\results'

colors = ['lightgreen', 'navy', 'cornflowerblue', 'red', 'lightcoral', 'firebrick']

def binom_ci_precision(proportion, nobs, method='beta', alpha=0.05):
    """
    Get precision for binomial proportion confidence interval
    """
    count = proportion * nobs
    ci = proportion_confint(count, nobs, method=method, alpha=alpha)
    ci_precision = ci[1] - proportion
    return ci_precision

def make_bar(data, error_bars='se', mean_amp=True, legend=False):
    
    use_cols = ['Fmax', 'cluster_05', 'cluster_01', 'BH', 'BY', 'BKY']
    if mean_amp:
        use_cols.insert(0, 'mean_amp')
    
    #Get values for error bars
    power_data = data.loc[:, use_cols].to_numpy().T
    if error_bars.lower() == 'se':
        stderr = np.sqrt( (power_data*(1-power_data)) / 10000 )
    elif error_bars.lower() == 'ci':
        stderr = binom_ci_precision(power_data, 10000)
    elif error_bars is None:
        stderr = None
    else:
        raise ValueError('Incorrect input for error_bars')
    
    #Plot
    labels = ['Fmax', 'Cluster 0.05', 'Cluster 0.01', 'BH FDR', 'BY FDR', 'BKY FDR']
    if mean_amp:
        labels.insert(0, 'Mean Amplitude')
        colors.insert(0, 'black')
    data.plot.bar(x='time_window', y=use_cols, label=labels, color=colors,
                  fontsize=12, yerr=stderr, legend=legend)
    plt.xticks(rotation='horizontal')
    plt.xlabel('')
    if legend:
        plt.legend(loc=(1.04,0), prop={'size': 12})
    

def make_power_figures():
    
    #Get all results csv files
    results_files = [file for file in os.listdir(results_dir) if file.endswith('.csv')]
    
    for results_file in results_files:
        
        #Load data
        data = pd.read_csv(join(results_dir, results_file))
        
        if 'Power' in results_file and 'Familywise' in results_file:
            
            if 'FamilywisePower' in results_file:
                mean_amp = True
            else:
                mean_amp = False
        
            #Make file with legend
            if not os.path.isfile(join(results_dir, 'legend.tif')):
                make_bar(data[0:3], legend=True)
                img_file = join(results_dir, 'legend.tif')
                plt.savefig(img_file, bbox_inches='tight', dpi=600)
                plt.close()
                
            #Make figures
            make_bar(data[0:3], error_bars='CI', mean_amp=mean_amp)
            img_file = join(results_dir, '%s_N400.tif' % results_file.strip('.csv'))
            plt.savefig(img_file, bbox_inches='tight', dpi=600)
            plt.close()
            
            make_bar(data[3:6], error_bars='CI', mean_amp=mean_amp) 
            img_file = join(results_dir, '%s_P300.tif' % results_file.strip('.csv'))
            plt.savefig(img_file, bbox_inches='tight', dpi=600)
            plt.close()
            
            make_bar(data[6:9], error_bars='CI', mean_amp=mean_amp)
            img_file = join(results_dir, '%s_P1.tif' % results_file.strip('.csv'))
            plt.savefig(img_file, bbox_inches='tight', dpi=600)
            plt.close()
            
def make_null_figures():
    pass

#%% Make EW figures

ew_files = [file for file in os.listdir(results_dir) if file.endswith('0-1000.csv')]

ew_file = ew_files[0]

data = pd.read_csv(join(results_dir, ew_files[0]))

#%%

bplot = data.loc[:, 'Fmax':].boxplot(whis=[5, 95], showfliers=False, return_type='dict')

for key in bplot.keys():
    i = 0
    for item in bplot[key]:
        item.set_linewidth(4)
        item.set_color(colors[int(i)])
        if key in ['whiskers', 'caps']:
            i += 0.5
        else:
            i += 1
        
