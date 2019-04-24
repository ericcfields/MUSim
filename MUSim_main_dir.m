%Return main directory for MUSim
%
%Author: Eric Fields
%Version Date: 22 March 2019
%
%Copyright (c) 2019, Eric C. Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.
%This software is provided "as is" and any express or implied warranties are disclaimed. 

function main_dir = MUSim_main_dir()
    if isunix()
        main_dir = '/gsfs0/data/fields/MUSim';
    else
        main_dir = 'C:\Users\ecfne\Documents\Eric\Research\Stats Simulations\MUSim';
    end
end
