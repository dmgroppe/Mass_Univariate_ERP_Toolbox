%Part 1
[types type_count]=list_event_types('syn02-s253-clean.set');

% Add bin information to the set files
if 0 % Only needs to be done once and is already done with the files I have
    bin_info2EEG('syn02-s253-clean.set','synVSnonsyn_blf.txt','syn02-s253-clean.set');
    for sub=[5 7 8 10],
        fname=['syn' num2str(sub,'%.2d') '-s253-clean.set'];
        bin_info2EEG(fname,'synVSnonsyn_blf.txt',fname);
    end
    
    for sub=[2 5 7 8 10],
        fname=['syn' num2str(sub,'%.2d') '-s254-clean.set'];
        bin_info2EEG(fname,'synVSnonsyn_blf.txt',fname);
    end
end

%% Create GND from EEGLAB set files (syn02-s253-clean.set,	syn02-s254-clean.set)
GND=sets2GND('gui','bsln',[-100 0],'exclude_chans',{'LO1','IO1','SO1','LO2'});
load synVSnonsyn.GND -MAT
gui_erp(GND)
GND
headinfo(GND)
plot_wave(GND,[1 2],'include_chans',{'FPz','Fz','Cz','Pz','Oz'});
help sets2GND

%% Create GND from ERPLAB erp files (s1_with_locs.erp,	s2_with_locs.erp)
pause; %%
GND=erplab2GND('gui','exclude_chans',{'VEOG','HEOG','A2'});
gui_erp(GND)
headinfo(GND)