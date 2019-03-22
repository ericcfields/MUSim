%Return an empty GND struct
%
%Author: Eric Fields
%Version Date: 22 March 2019

function GND = empty_GND()

    GND = struct;
    GND.exp_desc = '';
    GND.filename =  '';
    GND.filepath = '';
    GND.saved = 'no';
    GND.grands =  [];
    GND.grands_stde = [];
    GND.grands_t = [];
    GND.sub_ct = [];
    GND.chanlocs = [];
    GND.bin_info = [];
    GND.condesc = {};
    GND.time_pts = [];
    GND.bsln_wind = [];
    GND.odelay = [];
    GND.srate = 512;
    GND.indiv_fnames = {};
    GND.indiv_subnames = {};
    GND.indiv_traits = [];
    GND.indiv_bin_ct = [];
    GND.indiv_bin_raw_ct = [];
    GND.indiv_erps = [];
    GND.indiv_art_ics = {};
    GND.cals = [];
    GND.history = [];
    GND.t_tests = [];
    
end