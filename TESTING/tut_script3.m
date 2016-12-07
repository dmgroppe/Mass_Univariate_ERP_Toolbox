%% Part 3
GRP=GNDs2GRP('gui','create_difs','yes');
GRP
headinfo(GRP)
GRP=bin_opGRP(GRP,'(A+B)/n','Female',1,'Male',1,'Tone Counting Task (mean of Female and Male)');
gui_erp(GRP)
GRP=tmaxGRP(GRP,1,'time_wind',[25 250],'exclude_chans',{'lhz','rhz','lle','rle','A2'});
sig_raster(GRP,1);
gui_erp(GRP,'t_test',1);
%help tmaxGRP
GRP=tfdrGRP(GRP,1,'method','by','time_wind',[25 250],'exclude_chans',{'lhz','rhz','lle','rle','A2'}); 
%help tfdrGRP
GRP=clustGRP(GRP,1,'time_wind',[25 250],'exclude_chans',{'lhz','rhz','lle','rle','A2'},'chan_hood',.61,'thresh_p',.05);
%help clustGRP
GRP=tmaxGRP(GRP,1,'time_wind',[60 90; 95 150],'mean_wind','yes');
sig_raster(GRP,3);
sig_topo(GRP,4,'units','uV');
sig_topo(GRP,4,'units','t');
GRP=tfdrGRP(GRP,1,'time_wind',[60 90; 95 150],'mean_wind','yes','method','by');
GRP=clustGRP(GRP,1,'time_wind',[95 150],'mean_wind','yes','chan_hood',.61,'thresh_p',.05);
