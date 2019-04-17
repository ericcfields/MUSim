# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 15:28:57 2019

@author: Eric
"""

import os
from os.path import join
import pandas as pd

def make_csv(results_file, measure_name, measure_text, output_file=False, main_dir=None):
    
    if main_dir is None:
        main_dir = os.getcwd()
    
    with open(results_file) as f_in:
        results_text = f_in.readlines()
    
    out_df = pd.DataFrame(columns=('effect', 'time_window', 'electrodes'))
    
    r = -1 #data frame row index
    
    for i in range(len(results_text)):
        
        line = results_text[i]
        
        if 'SIMULATION SUMMARY' in line:
            r += 1
            
        elif 'Simulated experiments' in line:
            (field, value) = line.split('\t')
            out_df.loc[r, 'n_experiments'] = int(value)
            
        elif 'Permutations' in line:
            (field, value) = line.split('\t')
            out_df.loc[r, 'n_permutations'] = int(value)
            
        elif 'Sample size' in line:
            (field, value) = line.split('\t')
            out_df.loc[r, 'n_subjects'] = int(value)
            
        elif 'Effect:' in line:
            (field, value) = line.split('\t')
            out_df.loc[r, 'effect'] = os.path.basename(value).strip()
            
        elif 'Time window' in line:
            (field, value) = line.split('\t')
            out_df.loc[r, 'time_window'] = value.strip()
            
        elif 'Electrodes:' in line:
            (field, value) = line.split('\t')
            out_df.loc[r, 'electrodes'] = value.strip()
            
        elif 'Trials =' in line:
            (field, value) = line.split('\t')
            out_df.loc[r, 'n_trials'] = int(value.split()[0])
            
        elif 'MEAN WINDOW' in line:
            (field, value) = results_text[i+1].split('\t')
            out_df.loc[r, 'mean_amp'] = float(value)
            
        elif 'UNCORRECTED' in line:
            method = 'uncorrected'
        elif 'SIDAK' in line:
            method = 'sidak'
        elif 'FMAX' in line:
            method = 'Fmax'
        elif 'CLUSTER 0.05' in line:
            method = 'cluster_05'
        elif 'CLUSTER 0.01' in line:
            method = 'cluster_01'
        elif 'BH' in line:
            method = 'BH'
        elif 'BY' in line:
            method = 'BY'
        elif 'BKY' in line:
            method = 'BKY'
            
        elif measure_text in line:
            (field, value) = line.split('\t')
            out_df.loc[r, method] = float(value)
    
    if output_file:
        output_file = join(main_dir, 'results', 'MUSim_%s.csv' % measure_name)
        out_df.to_csv(output_file, index=False)
    
    return out_df

def main():
    main_dir = r'C:\Users\ecfne\Documents\Eric\Research\Stats Simulations\MUSim'
    results_file = join(main_dir, 'results', 'MUSim_power_results.txt')
    
    measures = {'FamilywisePower': 'Family-wise rejection rate across time points with effect (familywise power)',
                'FamilywsieTypeI': 'Family-wise rejection rate across time points with null effect (familywise Type I error)',
                'Mean_EW_Power': 'Mean rejection rate at individual time points with effect (element-wise power)',
                'Median_EW_Power': 'Median rejection rate at individual time points with effect (element-wise power)',
                'Mean_EW_TypeI': 'Mean rejection rate at individual time points with null effect (element-wise Type I error)',
                'Median_EW_TypeI': 'Median rejection rate at individual time points with null effect (element-wise Type I error)',
                'Mean_EW_FDR': 'Mean element-wise false discovery rate',
                'Median_EW_FDR': 'Median element-wise false discovery rate'}
    
    for measure_name in measures:
        make_csv(results_file, measure_name, measures[measure_name], output_file=True, main_dir=main_dir)
        
if __name__ == '__main__':
    main()
