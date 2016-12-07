%% Part2
load Visodbl_original.GND -MAT
GND.t_tests=[]; %not part of original script but sometimes Visodbl.GND has tests saved with it
headinfo(GND);
plot_wave(GND,[1 2]);
pause; %%
GND=bin_dif(GND,2,1,'X/O Task Target-Standards');
gui_erp(GND,'bin',5);
GND=tmaxGND(GND,5,'time_wind',[100 900],'output_file','visodbl_XO.txt');
pause; %%
%print -f1 -depsc odblXO_raster
sig_raster(GND,1);
gui_erp(GND,'t_test',1);
sig_wave(GND,1,'use_color','no');
%help tmaxGND
pause; %%
GND=tfdrGND(GND,5,'method','by','time_wind',[100 900],'output_file','p3xo_fdr.txt');
%help tfdrGND
pause; %%
GND=clustGND(GND,5,'time_wind',[100 900],'chan_hood',.61,'thresh_p',.05);
%help clustGND
pause; %%
GND=tmaxGND(GND,5,'time_wind',[150 250; 300 500],'mean_wind','yes');
sig_raster(GND,3);
sig_topo(GND,4,'units','uV');
sig_topo(GND,4,'units','t');
pause; %%
GND=tfdrGND(GND,5,'time_wind',[150 250; 300 500],'mean_wind','yes','method','by');
pause;
GND=clustGND(GND,5,'time_wind',[300 500],'chan_hood',.61,'thresh_p',.05,'mean_wind','yes');