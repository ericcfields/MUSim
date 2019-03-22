%Create template GND for MUSim
%
%Author: Eric Fields
%Version Date: 22 March 2019

eeglab;
clearvars; close all;

main_dir = MUSim_main_dir();

load(fullfile(main_dir, 'data', 'EmProb_13subs_128Hz.GND'), '-mat');

GND.exp_desc = 'SimulatedExp';
GND.filename = '';
GND.filepath = '';
GND.saved = 'no';
GND.grands = [];
GND.grands_stder = [];
GND.grands_t = [];
GND.sub_ct = [];
GND.chanlocs = GND.chanlocs(1:32);
GND.bin_info = [];
GND.bsln_wind = [];
GND.indiv_fnames = {};
GND.indiv_subnames = {};
GND.indiv_bin_ct = [];
GND.indiv_bin_raw_ct = [];
GND.indiv_erps = [];
GND.indiv_art_ics = {};
GND.history = {};
GND.F_tests = [];

GND = save_matmk(GND, 'MUSim_template.GND', fullfile(main_dir, 'data'), 1);


