%Return main directory for MUSim
%
%Author: Eric Fields
%Version Date: 22 March 2019

function main_dir = MUSim_main_dir()
	if isunix()
		main_dir = '/gsfs0/data/fields/MUSim';
    else
        main_dir = 'C:\Users\ecfne\Documents\Eric\Research\Stats Simulations\MUSim';
	end
end
