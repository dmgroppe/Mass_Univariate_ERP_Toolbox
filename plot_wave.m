% plot_wave() - Plots the ERPs from one or more bins in a GND or GRP 
%                variable in uV or t-scores.  Each channel is visualized 
%                on its own axis.
%
% Usage:
%  >> plot_wave(GND_GRP_or_fname,bins,varargin)
%
% Required Inputs:
%  GND_GRP_or_fname - A GND/GRP structure variable or the filename of
%                     a GND/GRP structure that has been saved to disk.
%                     To create a GND variable from Kutaslab ERP files (e.g.,
%                     *.nrm files) use avgs2GND.m.  To do the same from
%                     EEGLAB *.set files use sets2GND.m.  To create a
%                     a GRP structure use GNDs2GRP.m. See Mass Univariate ERP
%                     Toolbox documentation for detailed information about the
%                     format of GND and GRP variables. If you specifiy a filename
%                     be sure to include the file's path, unless the file is
%                     in the current working directory.
%  bins             - [vector] The index of one or more bins that you wish 
%                     to visualize. You can visualize up to eight bins.
%                     To see what test results are available, run 
%                     headinfo(GND) or headinfo(GRP).
%
% Optional Inputs:
%  exclude_chans - A cell array of channel labels to exclude from the
%                  visualization (e.g., {'A2','lle','rhe'}).  If you wish
%                  to exclude a single channel, you can enter it as a string 
%                  (e.g., 'A2') instead of a cell array. Use headinfo.m to see
%                  the channel labels stored in the GND/GRP variable. You cannot
%                  use both this option and 'include_chans' (below).{default: 
%                  not used, all channels visualized}
%  include_chans - A cell array of channel labels to include from the
%                  visualization (e.g., {'Fz','Cz','Pz'}).  If you wish
%                  to include a single channel, you can enter it as a string 
%                  (e.g., 'Cz') instead of a cell array. Use headinfo.m to see
%                  the channel labels stored in the GND/GRP variable. You cannot
%                  use both this option and 'exclude_chans' (above).{default: 
%                  not used, all channels visualized}
%  time_limits   - [min max] Time limits for the x-axes in units of ms
%                  (e.g., [-100 500]. {default: first and last time
%                  point}
%  time_ticks    - [vector] The time points (in units of ms) you want
%                  marked with ticks on the x-axes. {default: MATLAB
%                  default ticks}
%  tick_labels   - ['yes' or 'no'] If 'yes', major time tick marks will
%                  be labeled with their value (in milliseconds).
%                  {default: 'yes'}
%  cal_amp       - The amplitude (in microvolts or t-scores) of the 
%                  calibration bar drawn at time=0. {default: a third of
%                  waveform absolute maximum across all shown channels}
%  cal_width     - The width (in units of ms) of the "wings" on the top and
%                  bottom of the calibration bar. {default: scaled to time
%                  window}
%  wave_limits   - [min max] Limits for the y-axis in units of t-scores
%                  or uV (e.g., [-10 15]. {default: slightly bigger than
%                  min and max of waveform}
%  units         - ['uV' or 't'] If 'uV' topographies will be visualized
%                  in units of microvolts.  If 't', they will be shown
%                  in t-scores. {default: 'uV'}
%  legend        - ['none','corners' or 'box'] If 'none', no legend will be
%                  drawn. If 'corners', the descriptors for each bin will 
%                  be written in the corners of the figure.  The color of 
%                  the text will correspond to the waveforms it describes.
%                  If 'box' a standard MATLAB legend box will be drawn on 
%                  the figure indicating the bin descriptors for each 
%                  waveform. Note, if you use 'box,' you can manually move 
%                  the position of the legend by clicking-and-dragging the 
%                  legend box. {default: 'corners'}
%  ydir          - [1 or -1] If 1, positive is plotted up.  Otherwise,
%                  positive is represented as down. {default: -1}
%  fig_id        - [integer] The index # of the MATLAB figure in which
%                  the diagram will be produced.  Useful for overwriting
%                  old figures. {default: lowest unused index}
%  label_size    - [integer] The font size of the channel labels on top of the
%                  calibration bars.  If 0, channels labels are not drawn.
%                  {default: 8}
%  title         - [string] The plot title. {default: none}
%  verblevel     - An integer specifiying the amount of information you want
%                  this function to provide about what it is doing during runtime.
%                    Options are:
%                     0 - quiet, only show errors, warnings, and EEGLAB reports
%                     1 - stuff anyone should probably know
%                     2 - stuff you should know the first time you start working
%                         with a data set {default value}
%                     3 - stuff that might help you debug (show all
%                         reports)
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
% ?/??/2010-

%%%%%%%%%%%%%%%% FUTURE WORK %%%%%%%%%%%%%%%%%
%

function plot_wave(GND_GRP_or_fname,bins,varargin)

p=inputParser;
p.addRequired('GND_GRP_or_fname',@(x) isstruct(x) || ischar(x)); 
p.addRequired('bins',@isnumeric);
p.addParamValue('exclude_chans',[],@(x) ischar(x) || iscell(x));
p.addParamValue('include_chans',[],@(x) ischar(x) || iscell(x));
p.addParamValue('time_ticks',[],@isnumeric);
p.addParamValue('cal_amp',[],@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('cal_width',[],@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('fig_id',[],@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('verblevel',[],@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('time_limits',[],@(x) isnumeric(x) && (length(x)==2));
p.addParamValue('wave_limits',[],@(x) isnumeric(x) && (length(x)==2));
p.addParamValue('ydir',-1,@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('title',[],@(x) ischar(x));
p.addParamValue('units','uV',@ischar);
p.addParamValue('legend','corners',@(x) ischar(x) && (strcmpi('corners',x) || strcmpi('box',x) || strcmpi('none',x)));
p.addParamValue('tick_labels','yes',@ischar);
p.addParamValue('label_size',8,@(x) isnumeric(x) && (length(x)==1));

p.parse(GND_GRP_or_fname,bins,varargin{:});


% Manage VERBLEVEL
if isempty(p.Results.verblevel),
    VERBLEVEL=2; %not global, just local
else
    VERBLEVEL=p.Results.verblevel;
end

tick_labels=str2bool(p.Results.tick_labels);

%load GND/GRP variable if a filename was passed
if ischar(p.Results.GND_GRP_or_fname),
    if VERBLEVEL,
        fprintf('Loading GND struct from file %s.\n',p.Results.GND_GRP_or_fname);
    end
    load(p.Results.GND_GRP_or_fname,'-MAT');
    if ~exist('GND','var') && ~exist('GRP','var')
        error('File %s does not contain a GND or GRP variable.',p.Results.GND_GRP_or_fname);
    end
    if exist('GRP','var'),
        GND=GRP; %for simplicity GRP variables are re-named GND
        clear GRP;
    end
else
    GND=p.Results.GND_GRP_or_fname; %for simplicity GRP variables are re-named GND
end

%% Figure out which channels to ignore if any
%But first make sure exclude & include options were not both used.
if ~isempty(p.Results.include_chans) && ~isempty(p.Results.exclude_chans)
    error('You cannot use BOTH ''include_chans'' and ''exclude_chans'' options.');
end
if ischar(p.Results.exclude_chans),
    exclude_chans{1}=p.Results.exclude_chans;
elseif isempty(p.Results.exclude_chans)
    exclude_chans=[];
else
    exclude_chans=p.Results.exclude_chans;
end
if ischar(p.Results.include_chans),
    include_chans{1}=p.Results.include_chans;
elseif isempty(p.Results.include_chans)
    include_chans=[];
else
    include_chans=p.Results.include_chans;
end

n_chan=length(GND.chanlocs);
% exclude and include chan options
if ~isempty(exclude_chans),
    ignore_chans=zeros(1,length(exclude_chans)); %preallocate mem
    ct=0;
    for x=1:length(exclude_chans),
        found=0;
        for c=1:n_chan,
            if strcmpi(exclude_chans{x},GND.chanlocs(c).labels),
                found=1;
                ct=ct+1;
                ignore_chans(ct)=c;
            end
        end
        if ~found,
            watchit(sprintf('I attempted to exclude %s.  However no such electrode was found in GND variable.', ...
                exclude_chans{x}));
        end
    end
    ignore_chans=ignore_chans(1:ct);
    use_chans=setdiff(1:n_chan,ignore_chans);
elseif ~isempty(include_chans),
    use_chans=zeros(1,length(include_chans)); %preallocate mem
    ct=0;
    for x=1:length(include_chans),
        found=0;
        for c=1:n_chan,
            if strcmpi(include_chans{x},GND.chanlocs(c).labels),
                found=1;
                ct=ct+1;
                use_chans(ct)=c;
            end
        end
        if ~found,
            watchit(sprintf('I attempted to include %s.  However no such electrode was found in GND variable.', ...
                include_chans{x}));
        end
    end
    use_chans=use_chans(1:ct);
else
    use_chans=1:n_chan;
end

n_use_bins=length(bins);
if max(bins)>length(GND.bin_info)
    error('There are only %d bins in this GND/GRP variable, but you asked to plot Bin %d.\n', ...
        length(GND.bin_info),max(bins));
elseif min(bins)<1,
   error('All elements of the "bins" argument must be greater than 0.');  
end
if isempty(p.Results.fig_id)
    fig_h=figure;
else
    fig_h=figure(p.Results.fig_id); clf;
end
if n_use_bins>1,
    fig_title=[GND.exp_desc ' (Bins ' int2str(bins) ')'];
else
    fig_title=[GND.exp_desc ' (Bin ' int2str(bins) ')'];
end
set(fig_h,'name',fig_title,'paperpositionmode','auto');
%setting paperpositionmode to 'auto' means that if the figure is manually
%resized, the printed version of the figure will reflect the whatever the
%shown size was (at the time of printing)

%Waveform time limits
if isempty(p.Results.time_limits),
    time_lim=[GND.time_pts(1) GND.time_pts(end)];
else
    time_lim=p.Results.time_limits;
end

%Waveform amplitude limits
tpt(1)=find_tpt(time_lim(1),GND.time_pts);
tpt(2)=find_tpt(time_lim(2),GND.time_pts);
if strcmpi(p.Results.units,'uV'),
    dat=GND.grands(use_chans,tpt(1):tpt(2),bins);
else
    dat=GND.grands_t(use_chans,tpt(1):tpt(2),bins);
end
mx_wav=max(max(max(dat)));
mn_wav=min(min(min(dat)));
mx_abs=max(abs([mx_wav mn_wav]));
if isempty(p.Results.cal_amp),
    cal_amp=round(mx_abs/3);
else
    cal_amp=abs(p.Results.cal_amp);
end
if isempty(p.Results.wave_limits),
    y_lim(1)=min([mn_wav -cal_amp])*1.1;
    y_lim(2)=max([mx_wav cal_amp])*1.1;
else
    y_lim=p.Results.wave_limits;
end


%% Plot Waveforms
n_use_chans=length(use_chans);
n_tall=floor(sqrt(n_use_chans+1)); %add one axis for legend
n_wide=ceil((n_use_chans+1)/n_tall); %add one axis for legend
if n_use_bins>8,
   error('You are a crazy person! You can''t plot more than eight waveforms on a single plot.'); 
end
colors={'r','b','m','c','g','w','y'};
lwidth_plt=1;
if isempty(p.Results.cal_width)
    cal_width=(time_lim(2)-time_lim(1))/20;
else
    cal_width=p.Results.cal_width;
end
if isempty(p.Results.time_ticks)
    rng=time_lim(2)-time_lim(1);
    delt=rnd_orderofmag(rng/5);
    ticks=time_lim(1):delt:0;
    ticks=[ticks [0:delt:time_lim(2)]];
else
    ticks=p.Results.time_ticks;
end
tck_ht=cal_amp/5;
bin_ct=0;
bin_desc=cell(1,n_use_bins);
for b=bins,
    bin_ct=bin_ct+1;
    bin_desc{bin_ct}=GND.bin_info(b).bindesc;
    for c=1:n_use_chans,
        subplot(n_tall,n_wide,c);
        if b==bins(1),
            %first wave for this channel
            
            %wave amp=0 line
            h=plot([time_lim(1) time_lim(2)],[0 0],'k');
            set(h,'linewidth',lwidth_plt);
            hold on;
            %verticle part of cal bar
            h=plot([0 0],[-1 1]*cal_amp,'k');
            set(h,'linewidth',lwidth_plt);
            
            %horizontal parts of cal bar
            h=plot([-cal_width cal_width],[1 1]*cal_amp,'k');
            set(h,'linewidth',lwidth_plt);
            h=plot([-cal_width cal_width],[-1 -1]*cal_amp,'k');
            set(h,'linewidth',lwidth_plt);
   
            %time ticks
            tk_ct=0;
            for a=ticks,
                tk_ct=tk_ct+1;
                if mod(tk_ct,2)
                    %plot every other tick as half height
                    h=plot([a a],[-tck_ht tck_ht]*.5,'k');
                else
                    h=plot([a a],[-tck_ht tck_ht],'k');
                    last_lab=a; %for plotting time tick labels below
                end
                set(h,'linewidth',lwidth_plt);
            end
            
            %Channel label
            if p.Results.label_size>0,
                tt=text(0,sign(p.Results.ydir)*cal_amp*1.4,GND.chanlocs(use_chans(c)).labels);
                set(tt,'horizontalalignment','center','fontsize',p.Results.label_size, ...
                    'fontweight','bold');
            end
            
            set(gca,'visible','off');
            if p.Results.ydir<0
                set(gca,'ydir','reverse');
            end
        end
        
        if strcmpi(p.Results.units,'uV'),
            dat=GND.grands(use_chans(c),:,b);
        else
            dat=GND.grands_t(use_chans(c),:,b);
        end
        plot(GND.time_pts,dat,colors{bin_ct});
        
        if b==bins(end),
            %set limits
            axis([time_lim y_lim]);
        end
    end
    
    if strcmpi(p.Results.legend,'box') && (b==bins(end))
       [dummy obj_hndls]=legend('best',bin_desc);
       for bb=1:n_use_bins,
          set(obj_hndls(bb*2+n_use_bins-1),'color',colors{bb});
       end
    end
end

%% Plot empty axis for denoting tick marks and voltage size
subplot(n_tall,n_wide,n_use_chans+1);

if p.Results.label_size>0,
    unit_size=p.Results.label_size;
else
    unit_size=8;
end

%wave amp=0 line
h=plot([time_lim(1) time_lim(2)],[0 0],'k');
set(h,'linewidth',lwidth_plt);
hold on;
%verticle part of cal bar
h=plot([0 0],[-1 1]*cal_amp,'k');
set(h,'linewidth',lwidth_plt);

%horizontal parts of cal bar
h=plot([-cal_width cal_width],[1 1]*cal_amp,'k');
set(h,'linewidth',lwidth_plt);
h=plot([-cal_width cal_width],[-1 -1]*cal_amp,'k');
set(h,'linewidth',lwidth_plt);

%time ticks
tk_ct=0;
for a=ticks,
    tk_ct=tk_ct+1;
    if mod(tk_ct,2)
        %plot every other tick as half height
        h=plot([a a],[-tck_ht tck_ht]*.5,'k');
    else
        h=plot([a a],[-tck_ht tck_ht],'k');
    end
    set(h,'linewidth',lwidth_plt);
    
    %Tick labels
    if ~mod(tk_ct,2) && a,
        if tick_labels,
            if a==last_lab
                tt=text(a,-tck_ht*3,['     ' int2str(a) ' ms']);
            else
                tt=text(a,-tck_ht*3,int2str(a));
            end
            set(tt,'horizontalalignment','center','fontsize',unit_size, ...
                'fontweight','bold');
        end
    end
end

if p.Results.ydir<0
    set(gca,'ydir','reverse');
    sign_label='-';
else
    sign_label=[];
end

%Amplitude label
if strcmpi(p.Results.units,'uV'),
    unit_label='\muV';
else
    unit_label='t';
end
tt=text(0,sign(p.Results.ydir)*cal_amp*1.4,[sign_label num2str(cal_amp) ' ' unit_label]);
set(tt,'horizontalalignment','center','fontsize',unit_size,'fontweight','bold');

set(gca,'visible','off');
axis([time_lim y_lim]);

%Plot legend 
if strcmpi(p.Results.legend,'corners')
   leg_coorX=[.05 .05 .55 .55];
   leg_coorY=[.1 .05 .1 .05];
   bin_ct=0;
   for b=bins,
       bin_ct=bin_ct+1;
       tt=textsc(leg_coorX(bin_ct),leg_coorY(bin_ct),GND.bin_info(b).bindesc);
       set(tt,'horizontalalignment','left','fontsize',12,'color',colors{bin_ct});
   end
end

if ~isempty(p.Results.title),
   tt=textsc(p.Results.title,'title');
   set(tt,'fontweight','bold','fontsize',14);
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

