# -*- coding: utf-8 -*-
#Python 3.7, Pandas 0.24
"""
Compile MUSim simulation results into more usable csv files

Author: Eric Fields
Version Date: 17 April 2019
"""

import os
from os.path import join
import pandas as pd

def parse_results(results_file, measure_name, measure_text, output_file=False, main_dir=None):
    """
    Make csv for a particular measure
    
    INPUTS
    results_file  - Text file output by MUSim simulations
    measure_name  - Name of the measure for the output file
    measure_text  - Text of measure in results file
    output_file   - Boolean specifying whether to create csv output
    main_dir      - Main directory for MUSim
    
    OUTPUTS
    results_df   - data frame with parsed results
    """
    
    #Default to current directory as main directory
    if main_dir is None:
        main_dir = os.getcwd()
    
    #Read in results file
    with open(results_file) as f_in:
        results_text = f_in.readlines()
    
    #Set-up output data frame
    results_df = pd.DataFrame(columns=('effect', 'time_window', 'electrodes'))
    
    r = -1 #data frame row index
    
    #Parse lines and read relevant results into data frame
    for i in range(len(results_text)):
        
        line = results_text[i]
        
        #Starting to parse a new simulation
        if 'SIMULATION SUMMARY' in line:
            r += 1
        
        #Parse basic simulation parameters
        elif 'Simulated experiments' in line:
            (field, value) = line.split('\t')
            results_df.loc[r, 'n_experiments'] = int(value)
            
        elif 'Permutations' in line:
            (field, value) = line.split('\t')
            results_df.loc[r, 'n_permutations'] = int(value)
            
        elif 'Sample size' in line:
            (field, value) = line.split('\t')
            results_df.loc[r, 'n_subjects'] = int(value)
            
        elif 'Effect:' in line:
            (field, value) = line.split('\t')
            results_df.loc[r, 'effect'] = os.path.basename(value).strip()
            
        elif 'Time window' in line:
            (field, value) = line.split('\t')
            results_df.loc[r, 'time_window'] = value.strip()
            
        elif 'Electrodes:' in line:
            (field, value) = line.split('\t')
            results_df.loc[r, 'electrodes'] = value.strip()
            
        elif 'Trials =' in line:
            (field, value) = line.split('\t')
            results_df.loc[r, 'n_trials'] = int(value.split()[0])
            
        elif 'MEAN WINDOW' in line:
            (field, value) = results_text[i+1].split('\t')
            results_df.loc[r, 'mean_amp'] = float(value)
        
        #Find correction method for current point in file
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
        
        #Get the measure of interst and associate with method (from above)
        elif measure_text in line:
            (field, value) = line.split('\t')
            results_df.loc[r, method] = float(value)
    
    #Output results to csv if requested
    if output_file:
        output_file = join(main_dir, 'results', 'MUSim_%s.csv' % measure_name)
        results_df.to_csv(output_file, index=False)
    
    return results_df

def make_power_csvs():
    
    main_dir = r'C:\Users\ecfne\Documents\Eric\Research\Stats Simulations\MUSim'
    results_file = join(main_dir, 'results', 'MUSim_power_results.txt')
    
    #Measure names and texts
    measures = {'FamilywisePower': 'Family-wise rejection rate across time points with effect (familywise power)',
                'FamilywsieTypeI': 'Family-wise rejection rate across time points with null effect (familywise Type I error)',
                'Mean_EW_Power': 'Mean rejection rate at individual time points with effect (element-wise power)',
                'Median_EW_Power': 'Median rejection rate at individual time points with effect (element-wise power)',
                'Mean_EW_TypeI': 'Mean rejection rate at individual time points with null effect (element-wise Type I error)',
                'Median_EW_TypeI': 'Median rejection rate at individual time points with null effect (element-wise Type I error)',
                'Mean_EW_FDR': 'Mean element-wise false discovery rate',
                'Median_EW_FDR': 'Median element-wise false discovery rate'}
    
    #Produce csv for all measures
    for measure_name in measures:
        parse_results(results_file, measure_name, measures[measure_name], output_file=True, main_dir=main_dir)

def make_null_csvs():
    pass

def main():
    make_power_csvs()
    make_null_csvs()

if __name__ == '__main__':
    main()
