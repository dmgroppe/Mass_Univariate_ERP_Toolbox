% sig_wave()  - Illustrates results of a mass univariate t-test on a
%               butterfly plot of grand average ERPs, difference waves, or
%               spectra in voltage or t-scores.  Time windows used in the test 
%               are represented with vertical dashed lines. If plotting the
%               waves in units of t-scores, the critical t-score for
%               significance is represented with a horizontal dashed line.
%               Note that you can click on a wave to see the name of the
%               electrode where the wave was measured and the precise time
%               corresponding to that data point.  This information will
%               appear in a small box.  Click on the box again to make it
%               disappear.
%
% Usage:
%  >> sig_wave(GND_GRP_specGND_or_fname,test_id,varargin);
%
% Required Inputs:
%  GND_GRP_specGND_or_fname  - A GND/GRP/specGND structure variable or the 
%                      filename of such a variable that has been saved to
%                      disk. To create a GND variable from Kutaslab ERP 
%                      files (e.g., *.nrm files) use avgs2GND.m.  To do the
%                      same from EEGLAB *.set files use sets2GND.m.
%                      To create a GRP structure use GNDs2GRP.m. To create
%                      a specGND variable use sets2specGND.m. See Mass 
%                      Univariate ERP Toolbox documentation for detailed 
%                      information about the format of GND and GRP
%                      variables. If you specifiy a filename be sure to
%                      include the file's path, unless the file is in the
%                      current working directory.
%  test_id           - [integer] The index # of the t-test results
%                      stored in the GND/GRP/specGND variable that you wish
%                      to visualize. To see what test results are
%                      available, look at the t_tests field of the variable
%                      (e.g., GND.t_tests) or use headinfo.m or 
%                      headinfo_spec.m.
%
% Optional Inputs:
%  x_ticks          - [vector] The times/frequencies you want marked with ticks
%                     on the x-axis.  Note, because the EEG has been sampled 
%                     at discrete time points, not all times between the 
%                     minimum and maximum analyzed times/frequencies can be
%                     used.  When a tick is requested that is not available, 
%                     the closest possible value is used. If not
%                     specified, x_ticks are automatically chosen based
%                     on the time/frequency range covered by the diagram. 
%                     This option has no effect if the t-test was
%                     executed based on mean voltages/power within specified 
%                     time windows/frequency bands.
%  wave_ticks       - [vector] The waveform values you want marked with ticks
%                     on the y-axis. {default: MATLAB default ticks}
%  fig_id           - [integer] The index # of the MATLAB figure in which
%                     the diagram will be produced.  Useful for overwriting
%                     old figures. {default: lowest unused index}
%  units            - ['uV','dB', or 't'] If 'uV' topographies will be 
%                     visualized in units of microvolts.  If 'dB' spectral  
%                     power will be shown in decibels (10*log10((uV^2)/Hz)). 
%                     If 't', they will be shown in t-scores. {default:
%                     't'}
%  use_color        - [integer] If non-zero, non-gray scale colors will be
%                     used in the diagram. {default: 0}
%  x_limits         - [min max] Limits for the x-axis in units of ms for
%                     ERPs/difference waves (e.g., [-100 500]) and Hz for 
%                     power spectra (e.g., [9 13]). {default: first and 
%                     last time point or frequency}
%  wave_limits      - [min max] Limits for the y-axis in units of t-scores
%                     or uV (e.g., [-10 15]. {default: slightly bigger than
%                     min and max of waveform}
%  ydir             - [1 or -1] If 1, positive is plotted up for ERPs and
%                     difference waves.  Otherwise, positive is represented
%                     as down.  This argument has no effect if plotting 
%                     spectra. {default: -1}
%  test_line_width  - [integer] Width of the lines indicating test window
%                     boundaries and critical t-scores. {default: 2}
%  title            - [string] The axis title. {default: descriptor for the
%                     bin being visualized}
%  crit_val_box     - ['yes' or 'no'] If 'yes', a box indicating the alpha
%                     or q level for significance will be created next to
%                     the horizontal line indicating the critical t-score
%                     for significance.
%  verblevel        - An integer specifiying the amount of information you want
%                     this function to provide about what it is doing during runtime.
%                     Options are:
%                      0 - quiet, only show errors, warnings, and EEGLAB reports
%                      1 - stuff anyone should probably know
%                      2 - stuff you should know the first time you start working
%                          with a data set {default value}
%                      3 - stuff that might help you debug (show all
%                          reports)
%
% Global Variable:
%   VERBLEVEL - Mass Univariate ERP Toolbox level of verbosity (i.e., tells functions
%               how much to report about what they're doing during
%               runtime) set by the optional function argument 'verblevel'
%
% Notes:
% -The printed/exported figure will have the same dimensions as the figure
% on the display.  Thus you can undo figure clutter by re-sizing it.
%
%
% Author:
% David Groppe
% Kutaslab, 12/2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 12/15/2010-Compatible with specGND variables now too
%
% 6/4/2011-Throws error when user attemps to visualize cluster-based tests
%

%%%%%%%%%%%%%%%% FUTURE WORK %%%%%%%%%%%%%%%%%
%


function sig_wave(GND_GRP_specGND_or_fname,test_id,varargin)

p=inputParser;
p.addRequired('GND_GRP_specGND_or_fname',@(x) isstruct(x) || ischar(x));
p.addRequired('test_id',@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('x_ticks',[],@isnumeric);
p.addParamValue('wave_ticks',[],@isnumeric);
p.addParamValue('fig_id',[],@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('use_color',1,@(x) isnumeric(x) || ischar(x));
p.addParamValue('verblevel',[],@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('x_limits',[],@(x) isnumeric(x) && (length(x)==2));
p.addParamValue('wave_limits',[],@(x) isnumeric(x) && (length(x)==2));
p.addParamValue('ydir',-1,@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('test_line_width',2,@(x) isnumeric(x) && (length(x)==1) && x>0);
p.addParamValue('title',[],@(x) ischar(x));
p.addParamValue('units','t',@ischar);
p.addParamValue('crit_val_box','yes',@(x) ischar(x) && (strcmpi('yes',x) || strcmpi('no',x)));

p.parse(GND_GRP_specGND_or_fname,test_id,varargin{:});

% Manage VERBLEVEL
if isempty(p.Results.verblevel),
    VERBLEVEL=2; %not global, just local
else
    VERBLEVEL=p.Results.verblevel;
end

if ~(strcmpi(p.Results.units,'t') || strcmpi(p.Results.units,'uV') || strcmpi(p.Results.units,'dB'))
    error('Value of optional input ''units'' must be ''t'', ''uV'', or ''dB''.');
end

use_color=str2bool(p.Results.use_color);


%Load GND, GRP, or specGND struct
if ischar(GND_GRP_specGND_or_fname),
    fprintf('Loading GND, GRP, or specGND struct variable from file %s.\n',GND_GRP_specGND_or_fname);
    load(GND_GRP_specGND_or_fname,'-MAT');
    if ~exist('GND','var') && ~exist('GRP','var') && ~exist('specGND','var')
        error('File %s does not contain a GND, GRP, or specGND variable.',GND_GRP_specGND_or_fname);
    end
    if exist('GRP','var'),
        GND=GRP; %for simplicity GRP variables are re-named GND
        clear GRP;
    elseif exist('specGND','var'),
        GND=specGND; %for simplicity specGND variables are re-named GND
        clear specGND;
    end
    VerbReport(sprintf('Experiment: %s',GND.exp_desc),2,VERBLEVEL);
else
    GND=GND_GRP_specGND_or_fname; %for simplicity GRP and specGND variables are re-named GND
    clear GND_GRP_specGND_or_fname;
end

n_test=length(GND.t_tests);
if test_id>n_test,
   error('There are only %d t-tests stored with these data, but you requested Test #%d',n_test,test_id); 
end

%FDR or permutation test correction for multiple comparisons
if isnan(GND.t_tests(test_id).n_perm)
    fdr_crct=1;
    if VERBLEVEL,
        fprintf('Correcting for multiple comparisons via FDR procedure: %s\n', ...
            GND.t_tests(test_id).mult_comp_method);
    end
else
    fdr_crct=0;
    if VERBLEVEL,
        fprintf('Correcting for multiple comparisons via permutation test.\n');
    end
end

if isfield(GND.t_tests(1),'freq_band')
    %create temporary fields to make frequency domain plotting
    % easily compatible with time domain plot
    GND.t_tests(test_id).time_wind=GND.t_tests(test_id).freq_band;
    GND.t_tests(test_id).mean_wind=GND.t_tests(test_id).mean_band;
    GND.t_tests(test_id).used_tpt_ids=GND.t_tests(test_id).used_freq_ids;
    GND.grands_t=GND.grands_pow_dB_t;
    GND.grands=GND.grands_pow_dB;
    freq_step=GND.freqs(2)-GND.freqs(1);
    freq_ord=orderofmag(freq_step);
    ord_pow=log10(freq_ord); %useful for formatting x-tick labels
    if ord_pow>=0,
        n_dig_past_dot=0;
    else
        n_dig_past_dot=-ord_pow;
    end
    GND.time_pts=rnd_orderofmag(GND.freqs);
    freq_domain=1;
else
    freq_domain=0;
end

%Clean up code by copying some GND.t_tests variables to their own
%variable
use_bin=GND.t_tests(test_id).bin;
if strcmpi(GND.t_tests(test_id).mean_wind,'yes'),
    error('You can''t use sig_wave.m to visualize the results of t-tests performed on potentials averaged across a time window or frequency band.  Use sig_raster.m or sig_topo.m instead.');
elseif strcmpi(GND.t_tests(test_id).mult_comp_method,'cluster mass perm test'),
    error('You can''t use sig_wave.m to visualize the results of cluster based permutation tests.  Use sig_raster.m, gui_erp.m, or sig_topo.m instead.');
else
    if strcmpi(p.Results.units,'uV') || strcmpi(p.Results.units,'dB') ,
        grands_wav=squeeze(GND.grands(GND.t_tests(test_id).used_chan_ids, ...
            :,use_bin));
    else
        if GND.t_tests(test_id).null_mean,
            %recompute t-scores because mean of null hypothesis is not zero
            grands_wav=squeeze( (GND.grands(GND.t_tests(test_id).used_chan_ids, ...
                :,use_bin)-GND.t_tests(test_id).null_mean)./ ...
                GND.grands_stder(GND.t_tests(test_id).used_chan_ids, ...
                :,use_bin));
        else
            grands_wav=squeeze(GND.grands_t(GND.t_tests(test_id).used_chan_ids, ...
                :,use_bin));
        end
    end
end


if VERBLEVEL,
    if fdr_crct,
        fprintf('q level of critical t-scores: %f.\n',GND.t_tests(test_id).desired_alphaORq);
    else
        fprintf('Estimated alpha level of critical t-scores: %f.\n',GND.t_tests(test_id).estimated_alpha);
    end
end

if isempty(p.Results.fig_id)
    fig_h=figure;
else
    fig_h=figure(p.Results.fig_id); clf;
end
if freq_domain,
    set(fig_h,'name',['Bin ' int2str(use_bin) ' [' GND.bindesc{use_bin} ']'],'paperpositionmode','auto');
else
    set(fig_h,'name',['Bin ' int2str(use_bin) ' [' GND.bin_info(use_bin).bindesc ']'],'paperpositionmode','auto');
end
%setting paperpositionmode to 'auto' means that if the figure is manually
%resized, the printed version of the figure will reflect the whatever the
%shown size was (at the time of printing)


%% Plot waveforms
%First plot time=0 (if time domain), wave=0 lines
if isempty(p.Results.x_limits),
    time_lim=[GND.time_pts(1) GND.time_pts(end)];
else
    time_lim=p.Results.x_limits;
end
hd=plot(time_lim,[0 0],'k'); hold on; %time=0 line
set(hd,'linewidth',p.Results.test_line_width);
mx_abs_t=max(max(grands_wav));
mn_abs_t=min(min(grands_wav));
t_mx=mx_abs_t*1.1;
t_mn=mn_abs_t*1.1;
if isempty(p.Results.wave_limits),
    t_lim=[t_mn t_mx];
else
    t_lim=p.Results.wave_limits;
end
if freq_domain
    hd=plot([0 0],t_lim,'k'); % time=0 line
    set(hd,'linewidth',p.Results.test_line_width);
end


%plot waveforms
h_wav=plot(GND.time_pts,grands_wav');
axis([time_lim t_lim]);
chan_ct=0;
for w=h_wav',
    chan_ct=chan_ct+1;
    dat.lab=GND.chanlocs(GND.t_tests(test_id).used_chan_ids(chan_ct)).labels;
    dat.times=GND.time_pts;
    if freq_domain,
        bdf_code=['Cp = get(gca,''CurrentPoint''); ' ...
            'Xp=Cp(2,1);', ...
            'Yp=Cp(2,2);', ...
            'dat=get(gcbo,''userdata'');', ...
            'ht=text(Xp,Yp,[int2str(rnd_orderofmag(Xp)) '' Hz, '' dat.lab]);' ...
            'set(ht,''backgroundcolor'',''w'',''horizontalalignment'',''center'',''verticalalignment'',''middle'',''buttondownfcn'',''delete(gcbo);'',''edgecolor'',''k'');'];
    else
        bdf_code=['Cp = get(gca,''CurrentPoint''); ' ...
            'Xp=Cp(2,1);', ...
            'Yp=Cp(2,2);', ...
            'dat=get(gcbo,''userdata'');', ...
            'ht=text(Xp,Yp,[int2str(round(Xp)) '' ms, '' dat.lab]);' ...
            'set(ht,''backgroundcolor'',''w'',''horizontalalignment'',''center'',''verticalalignment'',''middle'',''buttondownfcn'',''delete(gcbo);'',''edgecolor'',''k'');'];
    end
    set(w,'ButtonDownFcn',bdf_code,'userdata',dat);
    if ~use_color,
        set(w,'color','k');
    end
end


%% plot test windows with critical vals
if use_color,
    crit_color='r';
else
    crit_color='k';
end

n_wind=size(GND.t_tests(test_id).time_wind,1);
for w=1:n_wind,
    hd=plot([1 1]*GND.t_tests(test_id).time_wind(w,1),t_lim,'k--');
    set(hd,'linewidth',p.Results.test_line_width);
    hd=plot([1 1]*GND.t_tests(test_id).time_wind(w,2),t_lim,'k--');
    set(hd,'linewidth',p.Results.test_line_width);
    if ~strcmpi(p.Results.units,'uV'),
        %don't plot critical values if waveform is in uV
        hd=plot(GND.t_tests(test_id).time_wind(w,:),[1 1]*GND.t_tests(test_id).crit_t(1),'r--');
        set(hd,'linewidth',p.Results.test_line_width,'color',crit_color);
        crit_t_for_text=GND.t_tests(test_id).crit_t(1); %t coordinate for writing alpha=.05 (or whatever)
        if length(GND.t_tests(test_id).crit_t)>1,
            hd=plot(GND.t_tests(test_id).time_wind(w,:),[1 1]*GND.t_tests(test_id).crit_t(2),'r--');
            set(hd,'linewidth',p.Results.test_line_width,'color',crit_color);
            if abs(t_lim(2))>abs(t_lim(1)),
                crit_t_for_text=GND.t_tests(test_id).crit_t(2);
            end
        end
    end
end

%% Axis labels
if strcmpi(p.Results.units,'t'),
    ht=ylabel('t-score');
else
    if freq_domain,
        ht=ylabel('10*log10(\muV^2/Hz)');
    else
        ht=ylabel('\muV');
        set(ht,'rotation',0);
    end
end
set(ht,'fontsize',14);
if freq_domain,
    ht=xlabel('Hz');
else
    ht=xlabel('Time (ms)');
end
set(ht,'fontsize',14);
if p.Results.ydir==-1 && ~freq_domain,
    set(gca,'ydir','reverse');
end

%% Title
if isempty(p.Results.title),
    if freq_domain,
        ht=title(GND.bindesc{use_bin});
    else
        ht=title(GND.bin_info(use_bin).bindesc);
    end
    set(ht,'fontsize',16);
else
    ht=title(p.Results.title);
    set(ht,'fontsize',16);
end

if strcmpi(p.Results.units,'t') && strcmpi(p.Results.crit_val_box,'yes'),
    %Add box denoting alpha level/FDR level
    
    %Figure out where to plot box
    tm_rng=time_lim(2)-time_lim(1);
    min_wind=min(min(GND.t_tests(test_id).time_wind));
    max_wind=max(max(GND.t_tests(test_id).time_wind));
    min_space=min_wind-time_lim(1);
    max_space=time_lim(2)-max_wind;
    if min_space>max_space,
        box_xloc=min_wind-tm_rng*.02;
        halign='right';
    else
        box_xloc=max_wind+tm_rng*.02;
        halign='left';
    end
    
    if strcmpi(GND.t_tests(test_id).mult_comp_method(1),'b')
        %FDR control used
        h_alph1=text(box_xloc,crit_t_for_text, ...
            ['q=' num2str(rnd_orderofmag(GND.t_tests(test_id).desired_alphaORq))]);
    else
        %FWER control used
        h_alph1=text(box_xloc,crit_t_for_text, ...
            ['\alpha=' num2str(rnd_orderofmag(GND.t_tests(test_id).estimated_alpha))]);
    end
    set(h_alph1,'color',crit_color,'fontweight','normal','fontsize',12, ...
        'horizontalalignment',halign,'backgroundcolor',[1 1 1], ...
        'edgecolor',[1 1 1]*.3,'clipping','on','fontname','fixedwidth');
end

%% Tick marks specified?
if ~isempty(p.Results.x_ticks),
    set(gca,'xtick',p.Results.x_ticks);
end
if ~isempty(p.Results.wave_ticks),
    set(gca,'ytick',p.Results.wave_ticks);
end


%
%% %%%%%%%%%%%%%%%%%%%%% function str2bool() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function bool=str2bool(str)
%function bool=str2bool(str)

if ischar(str),
    if strcmpi(str,'yes') || strcmpi(str,'y')
        bool=1;
    else
        bool=0;
    end
else
    bool=str;
end
