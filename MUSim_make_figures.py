# -*- coding: utf-8 -*-
"""
Make figures for MUSim paper

AUTHOR: Eric Fields
VERSION DATE: 26 June 2019
"""

import os
from os.path import join
import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportion_confint
import matplotlib.pyplot as plt


def binom_ci_precision(proportion, nobs, method='beta', alpha=0.05):
    """
    Get precision for binomial proportion confidence interval
    """
    count = proportion * nobs
    ci = proportion_confint(count, nobs, method=method, alpha=alpha)
    ci_precision = ci[1] - proportion
    return ci_precision

def make_power_bar(data, colors, error_bars='se', mean_amp=True, legend=False):
    
    use_colors = colors.copy()
    
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
    labels = ['Fmax', 'cluster (p≤0.05 threshold)', 'cluster (p≤0.05 threshold)', 
              'FDR (Benjamini & Hochberg, 1995)', 'FDR (Benjamini & Yekutieli, 2001)', 
              'FDR (Benjamini et al., 2006)']
    if mean_amp:
        labels.insert(0, 'mean amplitude')
        use_colors.insert(0, 'black')
    data.plot.bar(x='time_window', y=use_cols, label=labels, color=use_colors,
                  fontsize=16, yerr=stderr, legend=legend)
    plt.xticks(rotation='horizontal')
    plt.xlabel('')
    plt.ylim((0,1))
    if legend:
        plt.legend(loc=(1.04,0), prop={'size': 12})
    

def make_power_figures(colors, results_dir):
    
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
                make_power_bar(data[0:3], colors, legend=True)
                img_file = join(results_dir, 'legend.tif')
                plt.savefig(img_file, bbox_inches='tight', dpi=600)
                plt.close()
                
            #Make figures
            make_power_bar(data[0:3], colors, error_bars='CI', mean_amp=mean_amp)
            img_file = join(results_dir, '%s_N400.tif' % results_file.strip('.csv'))
            plt.savefig(img_file, bbox_inches='tight', dpi=600)
            plt.close()
            
            make_power_bar(data[3:6], colors, error_bars='CI', mean_amp=mean_amp) 
            img_file = join(results_dir, '%s_P300.tif' % results_file.strip('.csv'))
            plt.savefig(img_file, bbox_inches='tight', dpi=600)
            plt.close()
            
            make_power_bar(data[6:9], colors, error_bars='CI', mean_amp=mean_amp)
            img_file = join(results_dir, '%s_P1.tif' % results_file.strip('.csv'))
            plt.savefig(img_file, bbox_inches='tight', dpi=600)
            plt.close()
            
def make_null_figures(results_dir):

    #Get data
    data = pd.read_csv(join(results_dir, 'MUSim_Null_FamilywiseTypeI.csv'))
    data[['n_trials', 'n_subjects']] = data[['n_trials', 'n_subjects']].astype(int)
    
    #Plotting parameters
    use_cols = ['mean_amp', 'Fmax', 'cluster_05', 'cluster_01']
    labels = ['mean amplitude', 'Fmax', 'cluster (p ≤ 0.05 threshold)', 'cluster (p ≤ 0.01 threshold)']
    use_colors = ['black', 'lightgreen', 'navy', 'cornflowerblue']
    
    for time_wind in ('0 - 300', '300 - 1000'):
        for trials in (40, 20, 10):
            plot_subset = data[(data['time_window'] == time_wind) & (data['n_trials'] == trials)]
            proportions = plot_subset.loc[:, use_cols].to_numpy().T
            stderr = binom_ci_precision(proportions, 10000)
            
            #Make bar graph
            plot_subset.plot.bar(x='n_subjects', y=use_cols, label=labels, color=use_colors,
                          fontsize=16, yerr=stderr, legend=False)
            plt.xticks(rotation='horizontal')
            plt.xlabel('')
            plt.ylim((0,0.1))
            plt.axhline(y=0.05,linewidth=1, color='r', linestyle='--')
            plt.yticks(np.arange(1,11)/100)
            plt.xlabel('Number of Subjects', fontsize=18)
            
            #Save file
            img_file = join(results_dir, 'MUSim_Null_FamilywiseTypeI_%s_%dtrials.tif' % (time_wind, trials))
            plt.savefig(img_file, bbox_inches='tight', dpi=600)
            plt.close()
            
def make_EW_figures(colors, results_dir):

    ew_files = [file for file in os.listdir(results_dir) if 'Power_EW' in file and file.endswith('.csv')]
    
    for ew_file in ew_files:
        
        #Get data
        data = pd.read_csv(join(results_dir, ew_file))
        #Rename colums to labels to be used in figure
        data.columns = ['uncorrected', 'Sidak', 'Fmax', 'Clust0.05', 'Clust0.01', 'BH FDR', 'BY FDR', 'BKY FDR']
        
        #Make box plot
        bplot = data.loc[:, 'Fmax':].boxplot(whis=[5, 95], showfliers=False, 
                                             return_type='dict', patch_artist=True,
                                             fontsize=12)
        
        #For proporition measures, set standard y-scale
        if 'onset' not in ew_file and 'offset' not in ew_file:
            plt.ylim((0,1))
        
        #Update colors and line sizes
        for key in bplot.keys():
            i = 0
            for item in bplot[key]:
                item.set_linewidth(4)
                if key == 'medians':
                    item.set_color('black')
                else:
                    item.set_color(colors[int(i)])
                if key in ['whiskers', 'caps']:
                    i += 0.5
                else:
                    i += 1
        
        #Save figure
        img_file = join(results_dir, ew_file.strip('.csv') + '.tif')
        plt.savefig(img_file, bbox_inches='tight', dpi=600)
        plt.close()
        
def main():
    
    results_dir = r'C:\Users\ecfne\Documents\Eric\Research\Stats Simulations\MUSim\results'

    colors = ['lightgreen', 'navy', 'cornflowerblue', 'red', 'lightcoral', 'firebrick']
    
    make_power_figures(colors, results_dir)
    
    make_null_figures(results_dir)
    
    make_EW_figures(colors, results_dir)

if __name__ == '__main__':
    main()
