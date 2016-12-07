% sig_topo() - Plots scalp topographies of EEG efffects that have been
%              tested for significance in a family of tests.  Only works for
%              EEG effects derived from averaging time points within 
%              time windows (for tests of every single time point within
%              time windows use gui_erp.m or gui_pow.m).  Topographies can 
%              be visualized in microvolts (for ERPs), decibels (for spectral
%              power), or t-scores.
%             
% Usage:
%  >> sig_topo(GND_GRP_specGND_or_fname,test_id,varargin)
%
% Required Inputs:
%  GND_GRP_specGND_or_fname - A GND, GRP, or specGND structure variable or the 
%                             filename of such a variable that has
%                             been saved to disk. If you specifiy a filename be
%                             sure to include the file's path, unless the file is
%                             in the current working directory.
%  test_id          - [integer] The index # of the t-test results 
%                     stored in the GND/GRP/specGND variable that you wish to visualize.  
%                     To see what test results are available, look at 
%                     the "t_tests" field of your variable (e.g., GND.t_tests)
%                     or use the functions headinfo.m or headinfo_spec.m.
%
% Optional Inputs:
%  fig_id           - [integer] The index # of the MATLAB figure in which
%                     the diagram will be produced.  Useful for overwriting
%                     old figures. {default: lowest unused index}
%  units            - ['uV','dB', or 't'] If 'uV' topographies will be 
%                     visualized in units of microvolts.  If 'dB' spectral  
%                     power will be shown in decibels (10*log10((uV^2)/Hz)). 
%                     If 't', they will be shown in t-scores. {default: 't'}
%  title_on         - [0 or 1] If 1, the bin # and bin descriptor will be
%                     printed at the top of the figure. {default: 1}
%  one_scale        - [0 or 1] If 1, all topographies will be plot using
%                     the same color scale (useful for comparing effect 
%                     sizes across topographies).  If 0, each topography
%                     will be plot using a color scale tailored to it
%                     (provides maximal resolution of topography
%                     distribution).  {default: 1 if plotting in uV or dB, 
%                     0 if plotting t-scores}
%  scale_limits     - [min_value max_value] Minimum and maximum value of
%                     color scale used to illustrate topographies (e.g.,
%                     [-10 10]). If specified, one_scale parameter is set
%                     to 1. If not specified, color scales are set to +/-
%                     the maximum absolute value of the topography.
%                     {default: not used}                     
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
%   VERBLEVEL - Mass Univariate ERP Toolbox level of verbosity (i.e., tells 
%               functions how much to report about what they're doing during             
%               runtime) set by the optional function argument 'verblevel'
%
% Notes:
% -The printed/exported figure will have the same dimensions as the figure
% on the display.  Thus you can undo figure clutter by re-sizing it.
%
% -Note you need data at least three channels to interpolate a scalp topography.
%
% Author:
% David Groppe
% Kutaslab, 3/2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 5/6/2010-Made compatible with FDR control
%
% 12/14/2010-Compatible with specGND variables now too
%
% 9/10/2012-scale_limits option added

%%%%%%%%%%%%%%%% FUTURE WORK %%%%%%%%%%%%%%%%%
% -

function sig_topo(GND_GRP_specGND_or_fname,test_id,varargin)

p=inputParser;
p.addRequired('GND_GRP_specGND_or_fname',@(x) ischar(x) || isstruct(x));
p.addRequired('test_id',@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('units','t',@ischar);
p.addParamValue('fig_id',[],@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('title_on',1,@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('one_scale',[],@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('verblevel',[],@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('scale_limits',[],@(x) isnumeric(x) && (length(x)==2));
p.addParamValue('fontsize',12,@(x) isnumeric(x) && (length(x)==1));

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

if isfield(GND.t_tests(1),'freq_band')
    % We're dealing with spectral power analysis.  Create temporary fields 
    % to make frequency domain plotting compatible easily compatible with time domain plot
    GND.t_tests(test_id).time_wind=GND.t_tests(test_id).freq_band;
    GND.t_tests(test_id).mean_wind=GND.t_tests(test_id).mean_band;
    GND.t_tests(test_id).used_tpt_ids=GND.t_tests(test_id).used_freq_ids;
    GND.grands_t=GND.grands_pow_dB_t;
    GND.grands=GND.grands_pow_dB;
    freq_step=GND.freqs(2)-GND.freqs(1);
    freq_ord=orderofmag(freq_step);
    ord_pow=log10(freq_ord); %useful for formatting topo titles
    if ord_pow>=0,
        n_dig_past_dot=0;
    else
        n_dig_past_dot=-ord_pow;
    end
    GND.time_pts=rnd_orderofmag(GND.freqs);
    freq_domain=1;
    if strcmpi(p.Results.units,'uV'),
       watchit('You can''t plot spectral power topographies in units of microvolts.  They will be plot in decibels (i.e., 10*log10((uV^2)/Hz))) instead.');
    end
else
    freq_domain=0;
    if strcmpi(p.Results.units,'dB'),
        watchit('You can''t plot ERP topographies in units of decibels.  They will be plot in microvolts instead.');
    end
end

if ~strcmpi(GND.t_tests(test_id).mean_wind,'yes'),
   error('sig_topo.m only works for significance tests applied to data averaged across multiple time points.'); 
end

n_topos=length(GND.t_tests(test_id).used_tpt_ids);
wdth=floor(sqrt(n_topos));
ht=ceil(n_topos/wdth);

bin=GND.t_tests(test_id).bin;
chans=GND.t_tests(test_id).used_chan_ids;
if length(chans)<3,
   watchit('sig_topo.m cannot plot scalp topographies with less than 3 electrodes.');
   return
end
if isempty(p.Results.fig_id)
    fig_h=figure;
else
    fig_h=figure(p.Results.fig_id); clf;
end
if isfield(GND,'bin_info'),
    set(fig_h,'name',['Bin ' int2str(bin) ' [' GND.bin_info(bin).bindesc ']'],'paperpositionmode','auto');
else
    set(fig_h,'name',['Bin ' int2str(bin) ' [' GND.bindesc{bin} ']'],'paperpositionmode','auto');
end
%setting paperpositionmode to 'auto' means that if the figure is manually
%resized, the printed version of the figure will reflect whatever the
%shown size was (at the time of printing)

if strcmpi(p.Results.units,'uV') || strcmpi(p.Results.units,'dB') ,
    if freq_domain,
        units='dB';
    else
        units='\muV';
    end
    %Compute mean voltage or spectral power in each time window
    mns=zeros(length(chans),n_topos);
    for w=1:n_topos,
        mns(:,w)=mean(GND.grands(chans,GND.t_tests(test_id).used_tpt_ids{w},bin),2);
    end
    if ~isempty(p.Results.scale_limits)
        one_scale=1;
        maplimits=p.Results.scale_limits;
        VerbReport(sprintf('Using color scale limits of: %f to %f',maplimits(1),maplimits(2)),2,VERBLEVEL);
        VerbReport(sprintf('All topographies will have the same color scale.'),2,VERBLEVEL);
    elseif isempty(p.Results.one_scale) || (p.Results.one_scale<=0),
        one_scale=0;
        maplimits='absmax';
    else
        one_scale=1;
        mx=max(max(abs(mns)));
        maplimits=[-1 1]*mx;
    end
else
    units='t';
    if ~isempty(p.Results.scale_limits)
        one_scale=1;
        maplimits=p.Results.scale_limits;
        VerbReport(sprintf('Using color scale limits of: %f to %f',maplimits(1),maplimits(2)),2,VERBLEVEL);
        VerbReport(sprintf('All topographies will have the same color scale.'),2,VERBLEVEL);
    elseif isempty(p.Results.one_scale) || (p.Results.one_scale>0),
        one_scale=1;
        mx=max(max(abs(GND.t_tests(test_id).data_t)));
        maplimits=[-1 1]*mx;
    else
        one_scale=0;
        maplimits='absmax';
    end
end

if VERBLEVEL,
    if isnan(GND.t_tests(test_id).n_perm),
        fprintf('q level of critical t-scores: %f., FDR method: %s\n', ...
            GND.t_tests(test_id).desired_alphaORq,GND.t_tests(test_id).mult_comp_method);
    else
        fprintf('Estimated alpha level of critical t-scores: %f.\n', ...
            GND.t_tests(test_id).estimated_alpha);
    end
end

for a=1:n_topos,
    subplot(ht,wdth,a);
    time_wind=GND.t_tests(test_id).time_wind(a,:);
    if strcmpi(p.Results.units,'uV') || strcmpi(p.Results.units,'dB') ,
        mn=mns(:,a);
    else
        mn=GND.t_tests(test_id).data_t(:,a);
    end
    
    if isnan(GND.t_tests(test_id).n_perm),
        sig_chans=find(GND.t_tests(test_id).fdr_rej(:,a));
    else
        sig_chans=find(GND.t_tests(test_id).adj_pval(:,a)<GND.t_tests(test_id).estimated_alpha);
    end
    topoplotMK(mn,GND.chanlocs(chans),'emarker2',{sig_chans,'o',[1 1 1],4}, ...
        'maplimits',maplimits);
    if freq_domain,
        lab_form=['%.' int2str(n_dig_past_dot) 'f-%.' int2str(n_dig_past_dot)  'f Hz'];
    else
        lab_form='%d-%d ms';
    end
    
    hh=title(sprintf(lab_form,time_wind(1),time_wind(2)));
    set(hh,'fontsize',p.Results.fontsize);
    if ~one_scale,
        hcb=colorbar;
        set(gca,'fontsize',p.Results.fontsize);
        hy=ylabel(hcb,units);
        set(hy,'rotation',0,'verticalalignment','middle','position',[2 0 0], ...
            'fontsize',p.Results.fontsize);
        %set(hy,'rotation',0,'verticalalignment','middle','position',[4.9+n_topos*.2 0.03 1.0001]);
        %2.0310    0.0000         0
    end
end

if one_scale,
    hcb=cbar_nudge('vert',0,maplimits);
    hy=ylabel(hcb,units);
    set(hy,'rotation',0,'verticalalignment','middle');
end

% Plot title of entire figure
if p.Results.title_on,      
    if isfield(GND,'bin_info'),
        h=textsc(['Bin ' int2str(bin) ': ' GND.bin_info(bin).bindesc],'title');
    else
        h=textsc(['Bin ' int2str(bin) ': ' GND.bindesc{bin}],'title');
    end
    set(h,'fontsize',14,'fontweight','bold');
    if ~verLessThan('matlab','8')
        set(h,'position',[.5 .975 0]);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cbar_nudge (plots colorbar slightly to left of cbar.m
% search for "DG change" to find the line I edited
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handle]=cbar_nudge(arg,colors,minmax, grad)

if nargin < 2
  colors = 0;
end
posscale = 'off';
if nargin < 1
  arg = 'vert';
  ax = [];
else
  if isempty(arg)
    arg = 0;
  end
  if arg(1) == 0
    ax = [];
    arg = 'vert';
  elseif strcmpi(arg, 'pos')
    ax = [];
    arg = 'vert';
    posscale = 'on';
  else      
    if isstr(arg)
      ax = [];
    else
      ax = arg;
      arg = [];
    end
  end
end

if nargin>2
  if size(minmax,1) ~= 1 | size(minmax,2) ~= 2
    help cbar
    fprintf('cbar() : minmax arg must be [min,max]\n');
    return
  end
end
if nargin < 4
    grad = 5;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose colorbar position
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (length(colors) == 1) & (colors == 0)
  t = caxis;
else
  t = [0 1];
end
if ~isempty(arg)
  if strcmp(arg,'vert')  
    cax = gca;
    pos = get(cax,'Position');
    stripe = 0.04; 
    edge = 0.01;
    space = -.02;

    set(cax,'Position',[pos(1) pos(2) pos(3) pos(4)])
    rect = [pos(1)+pos(3)+space pos(2) stripe*pos(3) pos(4)];
    ax = axes('Position', rect);
  elseif strcmp(arg,'horiz')
    cax = gca;
    pos = get(cax,'Position');
    stripe = 0.075; 
    space = .1;  
    set(cax,'Position',...
        [pos(1) pos(2)+(stripe+space)*pos(4) pos(3) (1-stripe-space)*pos(4)])
    rect = [pos(1) pos(2) pos(3) stripe*pos(4)];
    ax = axes('Position', rect);
  end
else
  pos = get(ax,'Position');
  if pos(3) > pos(4)
    arg = 'horiz';
  else
    arg = 'vert';
  end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw colorbar using image()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

map = colormap;
n = size(map,1);

if length(colors) == 1
  if strcmp(arg,'vert')
      if strcmpi(posscale, 'on')
          image([0 1],[0 t(2)],[ceil(n/2):n-colors]');
      else
          image([0 1],t,[1:n-colors]');
      end;
      set(ax,'xticklabelmode','manual')
      set(ax,'xticklabel',[],'YAxisLocation','right')
      
  else
    image(t,[0 1],[1:n-colors]);
    set(ax,'yticklabelmode','manual')
    set(ax,'yticklabel',[],'YAxisLocation','right')
  end
  set(ax,'Ydir','normal','YAxisLocation','right')

else % length > 1

  if max(colors) > n
    error('Color vector excedes size of colormap')
  end
  if strcmp(arg,'vert')
    image([0 1],t,[colors]');
    set(ax,'xticklabelmode','manual')
    set(ax,'xticklabel',[])
  else
    image([0 1],t,[colors]);
    set(ax,'yticklabelmode','manual')
    set(ax,'yticklabel',[],'YAxisLocation','right')
  end  
  set(ax,'Ydir','normal','YAxisLocation','right')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust cbar ticklabels
%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 2 
  if strcmp(arg,'vert')
      Cax = get(ax,'Ylim');
  else
      Cax = get(ax,'Xlim');
  end;
  CBTicks = [Cax(1):(Cax(2)-Cax(1))/(grad-1):Cax(2)]; % caxis tick positions
  CBLabels = [minmax(1):(minmax(2)-minmax(1))/(grad-1):minmax(2)]; % tick labels
  
  dec = floor(log10(max(abs(minmax)))); % decade of largest abs value
  CBLabels = ([minmax]* [ linspace(1,0, grad);linspace(0, 1, grad)]);
  if dec<1
    CBLabels = round(CBLabels*10^(1-dec))*10^(dec-1);
  elseif dec == 1
    CBLabels = round(CBLabels*10^(2-dec))*10^(dec-2);
  else
    CBLabels = round(CBLabels);
  end

  if strcmp(arg,'vert')
    set(ax,'Ytick',CBTicks);
    set(ax,'Yticklabel',CBLabels);
  else
    set(ax,'Xtick',CBTicks);
    set(ax,'Xticklabel',CBLabels);
  end
end
handle = ax;

%%%%%%%%%%%%%%%%%%
% Adjust cbar tag
%%%%%%%%%%%%%%%%%%

set(ax,'tag','cbar')