# -*- coding: utf-8 -*-
#Python 3.7, Pandas 0.24
#
#Copyright (c) 2019, Eric C. Fields
#All rights reserved.
#This code is free and open source software made available under the 3-clause BSD license.
#This software is provided "as is" and any express or implied warranties are disclaimed. 
"""
Compile MUSim simulation results into more usable csv files

Author: Eric Fields
Version Date: 20 June 2019
"""

import os
from os.path import join
import pandas as pd

def parse_results(results_file, measure_text):
    """
    Read results into data frame
    
    INPUTS
    results_file  - Text file output by MUSim simulations
    measure_text  - Text of measure of interest in results file
    
    OUTPUTS
    results_df   - data frame with parsed results
    """
    
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
    
    return results_df

def make_power_csvs(results_file):
    
    #Measure names and texts
    measures = {'FamilywisePower': 'Family-wise rejection rate across time points with effect (familywise power)',
                'FamilywiseTypeI': 'Family-wise rejection rate across time points with null effect (familywise Type I error)',
                'FamilywiseFDR': 'Familywise FDR (proportion of sig studies that include false positive time point)',
                'TotalMissRate': 'Total miss rate (only null time points rejected)'}
    
    #Produce csv for all measures
    for measure_name in measures:
        results_df = parse_results(results_file, measures[measure_name])
        results_dir = os.path.dirname(results_file)
        output_file = join(results_dir, 'MUSim_Power_%s.csv' % measure_name)
        results_df.to_csv(output_file, index=False, float_format='%.3f')

def make_null_csvs(results_file):
    
    #Measure names and texts
    measures = {'FamilywiseTypeI': 'Family-wise rejection rate across time points with null effect (familywise Type I error)',
                'Mean_EW_TypeI': 'Mean rejection rate at individual time points with null effect (element-wise Type I error)',
                'Median_EW_TypeI': 'Median rejection rate at individual time points with null effect (element-wise Type I error)'}
    
    #Produce csv for all measures
    for measure_name in measures:
        results_df = parse_results(results_file, measures[measure_name])
        results_df.sort_values(['n_trials', 'n_subjects'], ascending=False, inplace=True)
        results_dir = os.path.dirname(results_file)
        output_file = join(results_dir, 'MUSim_Null_%s.csv' % measure_name)
        results_df.to_csv(output_file, index=False, float_format='%.3f')

def main():
    
    main_dir = r'C:\Users\ecfne\Documents\Eric\Research\Stats Simulations\MUSim'
    
    power_results_file = join(main_dir, 'results', 'MUSim_power_results.txt')
    make_power_csvs(power_results_file)
    
    null_results_file = join(main_dir, 'results', 'MUSim_null_results.txt')
    make_null_csvs(null_results_file)

if __name__ == '__main__':
    main()
