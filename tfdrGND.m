% tfdrGND() - Tests the null hypothesis that the grand average voltage
%             of a bin is mu or that the grand average within-subject 
%             difference between two bins is mu using one of several 
%             possible false discovery rate (FDR) procedures to control for 
%             mulitple comparisons (note, mu is assumed to be 0 by default).
%             This function requires individual subject ERPs to 
%             be stored in a "GND" structure and outputs the 
%             test results in a number of graphical and text formats.  
%             For analogous between-subject comparisons use the function 
%             tfdrGRP.m.
%             
% Usage:
%  >> [GND, p_values, data_t, crit_t, adj_p]=tfdrGND(GND_or_fname,bin,varargin)
%
% Required Inputs:
%   GND_or_fname - A GND structure variable or the filename of a 
%                  GND structure that has been saved to disk.  To 
%                  create a GND variable from Kutaslab ERP files (e.g.,
%                  *.mas files) use avgs2GND.m.  To do the same from 
%                  EEGLAB *.set files use sets2GND.m.  See Mass
%                  Univariate ERP Toolbox documentation for detailed 
%                  information about the format of a GND variable. If you 
%                  specifiy a filename be sure to include the file's path, 
%                  unless the file is in the current working directory.
%   bin          - [integer] The bin to contrast against a voltage of mu 
%                  Use the function headinfo.m to see what bins are stored 
%                  in a GND variable.  Use the function bin_dif.m to create
%                  a difference wave between two bins whose significance
%                  you can test with this function.
%
% Optional Inputs:
%   tail          - [1,0, or -1] An integer specifying the tail of the
%                   hypothesis test.  "1" means upper-tailed (i.e., alternative 
%                   hypothesis is that the ERP/difference wave is greater 
%                   than the mean of the null hypothesis).  "0" means two-
%                   tailed (i.e., alternative hypothesis is that the 
%                   ERP/difference wave is not equal to the mean of the null
%                   hypothesis). "-1" means lower-tailed (i.e., alternative 
%                   hypothesis is that the ERP/difference wave is less than
%                   the mean of the null hypothesis). {default: 0}
%   q             - A number between 0 and 1 specifying the family-wise
%                   q level of the test. q is the upper bound on the 
%                   expected proportion of rejected null hypotheses that are
%                   false rejections (i.e., the FDR). {default: 0.05}
%   method        - ['bh', 'by', or 'bky'] The procedure used to control
%                   the FDR. 'bh' is the classic Benjamini & Hochberg (1995)
%                   procedure, which is guaranteed to control FDR when the 
%                   tests are independent or positively dependent (e.g., 
%                   positively correlated Gaussians). 'by' is a much more
%                   conservative version of 'bh' that always controls FDR
%                   (regardless of the dependency structure of the tests--
%                   Benjamini & Yekutieli, 2001). 'bky' is a "two-stage"
%                   version of 'bh' that is more powerful than 'bh' when a 
%                   lot of the null hypotheses are false (Benjamini, Krieger, &
%                   Yekutieli, 2006).  'bky' is guaranteed to control FDR when the
%                   tests are independent and tends to be slightly less
%                   powerful than 'bh' when few or no null hypothese are
%                   false. {default: 'bh'}
%   time_wind     - 2D matrix of pairs of time values specifying the beginning 
%                   and end of one or more time windows in ms (e.g., 
%                   [160 180; 350 550]).  Every single time point in 
%                   the time window will be individually tested (i.e.,
%                   maximal temporal resolution) if mean_wind option is
%                   NOT used. Note, boundaries of time window(s) may not 
%                   exactly correspond to desired time window boundaries 
%                   because of temporal digitization (i.e., you only have
%                   samples every so many ms). {default: 0 ms to the end of
%                   the epoch}
%   time_block_dur- [integers] A number or numbers (in milliseconds) 
%                   specifying a duration of time blocks in which the time 
%                   windows specified by time_wind will be divided.  For
%                   example, if time_wind=[300 600] and time_block_dur=100,
%                   the 300 to 600 ms time window will be sub-divided into
%                   100 ms windows (i.e., time_wind will equal [300 396;
%                   400 496; 500 596]).  This is an easy way to break up
%                   larger time windows of interest into smaller windows
%                   for mean window analysis (if mean_wind option is set to
%                   'yes').  If you specify multiple time windows with
%                   time_wind, you can break them up using durations of
%                   different lengths.  For example, if time_wind=[150 250;
%                   400 900] and time_block_dur=[25 100], the first time
%                   window (150 to 250 ms) will be broken up into 25 ms
%                   windows and the second window (400 to 900 ms) will be
%                   broken up into 100 ms windows. {default: not used}
%   mean_wind     - ['yes' or 'no'] If 'yes', the t-tests will be performed 
%                   on the mean amplitude within each time window 
%                   specified by time_wind.  This sacrifices temporal 
%                   resolution to increase test power by reducing the number
%                   of comparisons.  If 'no', every single time point within
%                   time_wind's time windows will be tested individually.
%                   {default: 'no'}
%   null_mean     - [number] The mean of the null hypothesis (i.e., mu) in 
%                   units of microvolts. {default: 0}
%   exclude_chans - A cell array of channel labels to exclude from the
%                   t-tests (e.g., {'A2','lle','rhe'}).  This option 
%                   sacrifices spatial resolution to increase test power by 
%                   reducing the number of comparisons. Use headinfo.m to see
%                   the channel labels stored in the GND variable. You cannot
%                   use both this option and 'include_chans' (below).{default: 
%                   not used, all channels included in test}
%   include_chans - A cell array of channel labels to use in the t-tests
%                   (e.g., {'A2','lle','rhe'}).  All other channels will
%                   be ignored. This option sacrifices spatial resolution to 
%                   increase test power by reducing the number of comparisons.
%                   Use headinfo.m to see the channel labels stored in the GND
%                   variable. You cannot use both this option and 
%                   'exclude_chans' (above). {default: not used, all channels 
%                   included in test}
%   verblevel     - An integer specifiying the amount of information you want
%                   this function to provide about what it is doing during runtime.
%                    Options are:
%                      0 - quiet, only show errors, warnings, and EEGLAB reports
%                      1 - stuff anyone should probably know
%                      2 - stuff you should know the first time you start working
%                          with a data set {default value}
%                      3 - stuff that might help you debug (show all
%                          reports)
%   plot_gui      - ['yes' or 'no'] If 'yes', a GUI is created for
%                   visualizing the results of the t-tests using the 
%                   function gui_erp.m. The GUI vizualizes the grand average 
%                   ERPs in each bin via various stats (uV, t-scores), shows 
%                   topographies at individual time points, and illustrates 
%                   which electrodes significantly differ from the null 
%                   hypothesis.  This option does not work if mean_wind 
%                   option is set to 'yes.' This GUI can be reproduced using
%                   the function gui_erp.m. {default: 'yes'}
%   plot_raster   - ['yes' or 'no'] If 'yes', a two-dimensional (time x channel)
%                   binary "raster" diagram is created to illustrate the
%                   results of the t-tests.  Significant negative and
%                   positive deviations from the null hypothesis are shown
%                   as black and white rectangles respectively. Non-
%                   significant comparisons are shown as gray rectangles. 
%                   Clicking on the rectangles will show you the 
%                   corresponding time and channel label for that
%                   rectangle. This figure can be reproduced with the 
%                   function sig_raster.m. {default: 'yes'}
%   plot_mn_topo  - ['yes' or 'no'] If 'yes', the topographies of the mean
%                   voltages/effects in each time window are produced.  More
%                   specifically, two figures are produced: one showing the
%                   topographies in uV the other in t-scores. Significant/
%                   nonsignificant comparisons are shown as white/black 
%                   electrodes. Clicking on electrodes will show the
%                   electrode's name.  This figure can be reproduced with
%                   the function sig_topo.m.  This option has NO effect if
%                   mean_wind option is set to 'no'. {default: 'yes'}
%   output_file   - A string indicating the name of a space delimited text
%                   file to produce containing the p-values of all comparisons 
%                   and the details of the test (e.g., FDR method used, 
%                   family-wise q level, etc...). If mean_wind option is
%                   set to 'yes,' t-scores of each comparison are also 
%                   included since you cannot derive them from the t-scores
%                   at each time point/electrode in a simple way. When 
%                   importing this file into a spreadsheet be sure NOT to 
%                   count consecutive spaces as multiple delimiters. If bh
%                   or bh FDR control procedures are used, FDR adjusted 
%                   p-values (also called "q-values") will be output to the 
%                   text file.  If method bky is used, unadjusted p-values 
%                   will be output since it is not clear how to compute FDR
%                   adjusted p-values for this method. {default: none}
%   save_GRP      - ['yes' or 'no'] If 'yes', the GRP variable will be
%                   saved to disk after the t-tests have been completed 
%                   and added to it. User will first be prompted to verify 
%                   file name and path. {default: 'yes'}
%
% Outputs:
%   GND           - GND structure variable.  This is the same as
%                   the input GND variable with one addition: the 
%                   field GND.t_tests will contain the results of the 
%                   t-tests and the test parameters. 
%   p_values      - A two-dimensional matrix (channel x time) of the
%                   p-values of each comparison (no correction for multiple
%                   comparisons).
%   data_t        - A two-dimensional matrix (channel x time) of the
%                   t-scores of each comparison.
%   crit_t        - The critical t-score(s) for the test.  Any t-scores
%                   that are more extreme than the critical t-score(s) 
%                   significantly deviate from 0.
%   adj_p         - FDR corrected p-values (also called q-values). Note, 
%                   FDR corrected p-values can be greater than 1. For bky 
%                   procedure adj_p is NaN since it is not clear how to
%                   compute adjusted p-values for this procedure.
%                   
%
%   Note also that a great deal of information about the test is displayed 
%   in the MATLAB command window.  You can easiy record of all this
%   information to a text file using the MATLAB command "diary."
%
% Global Variables:
%   VERBLEVEL = Mass Univariate ERP Toolbox level of verbosity (i.e., tells 
%               functions how much to report about what they're doing during
%               runtime) set by the optional function argument 'verblevel'
%
% Notes:
% -To add a difference wave to a GND variable, use the function "bin_dif.m".
%
% References:
%   Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery 
% rate: A practical and powerful approach to multiple testing. Journal of 
% the Royal Statistical Society. Series B (Methodological), 57(1), 289-300.
%
%   Benjamini, Y., Krieger, A. M., & Yekutieli, D. (2006). Adaptive linear 
% step-up procedures that control the false discovery rate. Biometrika, 
% 93(3), 491-507. 
%
%   Benjamini, Y., & Yekutieli, D. (2001). The control of the false 
% discovery rate in multiple testing under dependency. The Annals of 
% Statistics, 29(4), 1165-1188. 
%
% Author:
% David Groppe
% Kutaslab, 5/2010


%%%%%%%%%%%%%%%% FUTURE WORK %%%%%%%%%%%%%%%%%
% 

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 4/13/2011-FDR corrected p-values added for by and bh procedures. I don't
% know if there's a way to compute them for bky procedure.


function [GND, p_values, data_t, crit_t, adj_p]=tfdrGND(GND_or_fname,bin,varargin)

global VERBLEVEL;

p=inputParser;
p.addRequired('GND_or_fname',@(x) ischar(x) || isstruct(x));
p.addRequired('bin',@(x) isnumeric(x) && (length(x)==1) && (x>0));
p.addParamValue('tail',0,@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('q',0.05,@(x) isnumeric(x) && (x>0) && (x<1));
p.addParamValue('method','bh',@(x) strcmpi(x,'bh') || strcmpi(x,'by') || strcmpi(x,'bky'));
p.addParamValue('time_wind',[],@(x) isnumeric(x) && (size(x,2)==2));
p.addParamValue('time_block_dur',[],@isnumeric);
p.addParamValue('mean_wind','no',@(x) ischar(x) && (strcmpi(x,'yes') || strcmpi(x,'no')));
p.addParamValue('null_mean',0,@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('exclude_chans',[],@(x) ischar(x) || iscell(x));
p.addParamValue('include_chans',[],@(x) ischar(x) || iscell(x));
p.addParamValue('verblevel',[],@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('plot_gui','yes',@(x) ischar(x) && (strcmpi(x,'yes') || strcmpi(x,'no')));
p.addParamValue('plot_raster','yes',@(x) ischar(x) && (strcmpi(x,'yes') || strcmpi(x,'no')));
p.addParamValue('plot_mn_topo',[],@(x) ischar(x) && (strcmpi(x,'yes') || strcmpi(x,'no')));
p.addParamValue('output_file',[],@ischar);
p.addParamValue('save_GND','yes',@(x) ischar(x) && (strcmpi(x,'yes') || strcmpi(x,'no')));

p.parse(GND_or_fname,bin,varargin{:});


if isempty(p.Results.verblevel),
    if isempty(VERBLEVEL),
        VERBLEVEL=2;
    end
else
   VERBLEVEL=p.Results.verblevel; 
end
    
mean_wind=str2bool(p.Results.mean_wind);

%Load GND struct
if ischar(GND_or_fname),
    fprintf('Loading GND struct from file %s.\n',GND_or_fname);
    load(GND_or_fname,'-MAT');
    if ~exist('GND')
        error('File %s does not contain a GND variable.',GND_or_fname);
    end
else
    GND=GND_or_fname;
    clear GND_or_fname;
    fldnames=fieldnames(GND);
    if ismember('group_desc',fldnames),
       error('You passed a GRP variable to this function instead of a GND variable.'); 
    end
end
[n_chan, n_pt, n_bin, total_subs]=size(GND.indiv_erps);
VerbReport(sprintf('Experiment: %s',GND.exp_desc),2,VERBLEVEL);

if (bin>n_bin),
    error('There is no Bin %d in this GND variable.',bin); 
end

%Use only subs with data in relevant bin(s)
use_subs=find(GND.indiv_bin_ct(:,p.Results.bin));
n_sub=length(use_subs);
VerbReport(sprintf('%d out of %d participants have data in relevant bin.',n_sub,total_subs), ...
    1,VERBLEVEL);


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


%% Find time points
if isempty(p.Results.time_wind),
    time_wind=[0 GND.time_pts(end)]; %default time window
else
    time_wind=p.Results.time_wind;
end
time_wind=sort(time_wind,2); %first make sure earlier of each pair of time points is first
time_wind=sort(time_wind,1); %next sort time windows from earliest to latest onset
n_wind=size(time_wind,1);
if ~isempty(p.Results.time_block_dur),
    tbd=p.Results.time_block_dur;
    n_block=length(tbd);
    if (n_block>1) && (n_block~=n_wind),
       error('If you specify more than one time block duration you need to provide exactly one duration for every time window (in this case %d).',n_wind); 
    elseif (n_block==1) && (n_wind>1),
        tbd=repmat(tbd,1,n_wind);
    end
        
    t_winds=[];
    one_step=1000/GND.srate;
    for a=1:n_wind,
        wind_strt=time_wind(a,1);
        wind_stop=time_wind(a,2);
        new_winds=[wind_strt:tbd(a):wind_stop];
        for b=1:(length(new_winds)-1),
            t_winds=[t_winds; new_winds(b) new_winds(b+1)-one_step];
        end
    end
    time_wind=t_winds;
    clear t_winds;
    n_wind=size(time_wind,1);
end

if mean_wind,
    use_tpts=cell(1,n_wind);
else
    use_tpts=[];
end
for a=1:n_wind,
    VerbReport(sprintf('Time Window #%d:',a),1,VERBLEVEL);
    VerbReport(sprintf('Attempting to use time boundaries of %d to %d ms for hypothesis test.',time_wind(a,1),time_wind(a,2)), ...
        1,VERBLEVEL);
    start_tpt=find_tpt(time_wind(a,1),GND.time_pts);
    end_tpt=find_tpt(time_wind(a,2),GND.time_pts);
    if mean_wind,
        use_tpts{a}=[start_tpt:end_tpt];
    else
        use_tpts=[use_tpts [start_tpt:end_tpt]];
    end
    %replace desired time points with closest matches
    time_wind(a,1)=GND.time_pts(start_tpt);
    time_wind(a,2)=GND.time_pts(end_tpt);
    VerbReport(sprintf('Exact window boundaries are %d to %d ms (that''s from time point %d to %d).', ...
        time_wind(a,1),time_wind(a,2),start_tpt,end_tpt),1,VERBLEVEL);
end
if ~mean_wind,
    use_tpts=unique(use_tpts); %sorts time points and gets rid of any redundant time points
end

%% Compile data
if mean_wind,
    %Take mean amplitude in time blocks and then test
    erps=zeros(length(use_chans),n_wind,n_sub);
    for a=1:n_wind,
        for sub=1:n_sub,
            erps(:,a,sub)=mean(GND.indiv_erps(use_chans,use_tpts{a},bin,use_subs(sub)),2);
        end
    end
else
    %Use every single time point in time window(s)
    n_use_tpts=length(use_tpts);
    erps=zeros(length(use_chans),n_use_tpts,n_sub);    
    for sub=1:n_sub,
        erps(:,:,sub)=GND.indiv_erps(use_chans,use_tpts,bin,use_subs(sub));
    end
end


%% Report tail of test & q levels
VerbReport(sprintf('Testing null hypothesis that the grand average ERPs in Bin %d (%s) have a mean of %f microvolts.',bin, ...
    GND.bin_info(bin).bindesc,p.Results.null_mean),1,VERBLEVEL);
if p.Results.tail==0
    VerbReport(sprintf('Alternative hypothesis is that the ERPs differ from %f (i.e., two-tailed test).',p.Results.null_mean), ...
        1,VERBLEVEL);
elseif p.Results.tail<0,
    VerbReport(sprintf('Alternative hypothesis is that the ERPs are less than %f (i.e., lower-tailed test).',p.Results.null_mean), ...
        1,VERBLEVEL);
else
    VerbReport(sprintf('Alternative hypothesis is that the ERPs are greater than %f (i.e., upper-tailed test).',p.Results.null_mean), ...
        1,VERBLEVEL);
end


[p_values, data_t]=fast_t1(erps-p.Results.null_mean,p.Results.tail,VERBLEVEL);
switch lower(p.Results.method)
    case 'bh'
        VerbReport('FDR control procedure: Benjamini & Hochberg (independent or positive dependency)',1,VERBLEVEL);
        [h_rej, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_values,p.Results.q,'pdep','no');
    case 'by'
        VerbReport('FDR control procedure: Benjamini & Yekutieli (arbitrary dependency)',1,VERBLEVEL);
        [h_rej, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_values,p.Results.q,'dep','no');
    case 'bky'
        VerbReport('FDR control procedure: Benjamini, Krieger, & Yekutieli (two-stage)',1,VERBLEVEL);
        [h_rej, crit_p]=fdr_bky(p_values,p.Results.q,'no');
        adj_p=NaN;
end

sig_tpts=find(sum(h_rej));
if isempty(sig_tpts),
    if p.Results.tail==0,
        crit_t=[NaN NaN];
    else
        crit_t=NaN;
    end    
    fprintf('ERPs are NOT significantly different from %f (q=%f) at any time point/window analyzed.\n', ...
        p.Results.null_mean,p.Results.q);
    if ~strcmpi(p.Results.method,'bky'),
        %don't have adjusted p-values for bky method
        fprintf('All FDR adjusted p-values>=%f\n',min(min(p_values)));
    end
else
    crit_t=min(abs(data_t(crit_p==p_values)));
    if p.Results.tail<0,
        crit_t=-crit_t;
    elseif p.Results.tail==0,
        crit_t=[-crit_t crit_t];
    end
    
    %Command line summary of results
    VerbReport(['Critical t-score(s):' num2str(crit_t)],1,VERBLEVEL);
    if p.Results.tail
        %one-tailed test
        tw_alpha=1-cdf('t',max(abs(crit_t)),n_sub-1);
    else
        %two-tailed test
        tw_alpha=(1-cdf('t',max(abs(crit_t)),n_sub-1))*2;
    end
    VerbReport(sprintf('That corresponds to a test-wise alpha level of %f.',tw_alpha),1,VERBLEVEL);
    VerbReport(sprintf('Bonferroni test-wise alpha would be %f.',p.Results.q/numel(p_values)) ...
        ,1,VERBLEVEL);
    
    fprintf('Total number of significant differences: %d\n',sum(sum(h_rej)));
    fprintf('Estimated upper bound on expected number of false discoveries: %.1f\n',sum(sum(h_rej))*p.Results.q);
    fprintf('Significant differences from zero (in order of earliest to latest):\n');
    
    if ~strcmpi(p.Results.method,'bky'),
        %don't have adjusted p-values for bky method
        max_sig_p=0;
        min_sig_p=max(max(adj_p));
    end
    for t=sig_tpts,
        if mean_wind
            %time windows instead of time points
            fprintf('%d to %d ms, electrode(s): ',GND.time_pts(use_tpts{t}(1)), ...
                GND.time_pts(use_tpts{t}(end)));
        else
            fprintf('%d ms, electrode(s): ',GND.time_pts(use_tpts(t)));
        end
        sig_elec=find(h_rej(:,t));
        ct=0;
        for c=sig_elec',
            ct=ct+1;
            if ~strcmpi(p.Results.method,'bky'),
                %don't have adjusted p-values for bky method
                if adj_p(c,t)>max_sig_p,
                    max_sig_p=adj_p(c,t);
                end
                if adj_p(c,t)<min_sig_p,
                    min_sig_p=adj_p(c,t);
                end
            end
            if ct==length(sig_elec),
                fprintf('%s.\n',GND.chanlocs(use_chans(c)).labels);
            else
                fprintf('%s, ',GND.chanlocs(use_chans(c)).labels);
            end
        end
    end
    if ~strcmpi(p.Results.method,'bky'),
        fprintf('All significant corrected p-values are between %f and %f\n',max_sig_p,min_sig_p);
    end
end


%Add t-test results to GND struct
n_t_tests=length(GND.t_tests);
neo_test=n_t_tests+1;
GND.t_tests(neo_test).bin=bin;
GND.t_tests(neo_test).time_wind=time_wind;
GND.t_tests(neo_test).used_tpt_ids=use_tpts;
n_use_chans=length(use_chans);
include_chans=cell(1,n_use_chans);
for a=1:n_use_chans,
   include_chans{a}=GND.chanlocs(use_chans(a)).labels; 
end
GND.t_tests(neo_test).include_chans=include_chans;
GND.t_tests(neo_test).used_chan_ids=use_chans;
GND.t_tests(neo_test).mult_comp_method=p.Results.method;
GND.t_tests(neo_test).n_perm=NaN;
GND.t_tests(neo_test).desired_alphaORq=p.Results.q;
GND.t_tests(neo_test).estimated_alpha=NaN;
GND.t_tests(neo_test).null_mean=p.Results.null_mean;
if mean_wind,
    GND.t_tests(neo_test).data_t=data_t;
    GND.t_tests(neo_test).mean_wind='yes';
else
    GND.t_tests(neo_test).data_t='See GND.grands_t';
    GND.t_tests(neo_test).mean_wind='no';
end
GND.t_tests(neo_test).crit_t=crit_t;
GND.t_tests(neo_test).df=length(use_subs)-1;
GND.t_tests(neo_test).adj_pval=adj_p; % is NaN for BKY FDR procedure
GND.t_tests(neo_test).fdr_rej=h_rej;
GND.t_tests(neo_test).seed_state=NaN;
GND.t_tests(neo_test).clust_info=NaN;
GND.t_tests(neo_test).chan_hood=NaN;

if strcmpi(p.Results.plot_raster,'yes'),
    sig_raster(GND,neo_test,'verblevel',0,'use_color','rgb');
end

if mean_wind,
    if strcmpi(p.Results.plot_mn_topo,'yes') || isempty(p.Results.plot_mn_topo),
        sig_topo(GND,neo_test,'units','t','verblevel',0); %t-score topographies
        sig_topo(GND,neo_test,'units','uV','verblevel',0); %microvolt topographies
    end
else
    if strcmpi(p.Results.plot_gui,'yes'),
        gui_erp('initialize','GNDorGRP',GND,'t_test',neo_test,'stat','t','verblevel',1);
    end
end

%% Write results of test to text file if requested
if ~isempty(p.Results.output_file)
    [fid msg]=fopen(p.Results.output_file,'w');
    if fid==-1,
        error('Cannot create file %s for writing.  According to fopen.m: %s.', ...
            p.Results.file,msg);
    else
        %Write header column of times
        % Leave first column blank for channel labels   
        if mean_wind,
            for t=1:n_wind
                fprintf(fid,' %d-%d',GND.time_pts(use_tpts{t}(1)), ...
                    GND.time_pts(use_tpts{t}(end)));
            end
            
            %write a couple spaces and then write header for t-scores
            fprintf(fid,'  ');
            for t=1:n_wind
                fprintf(fid,' %d-%d',GND.time_pts(use_tpts{t}(1)), ...
                    GND.time_pts(use_tpts{t}(end)));
            end
        else
            for t=use_tpts,
                fprintf(fid,' %d',GND.time_pts(t));
            end
        end
        fprintf(fid,' Milliseconds\n');
        
        % Write channel labels and p-values
        chan_ct=0;
        for c=use_chans,
            chan_ct=chan_ct+1;
            fprintf(fid,'%s',GND.chanlocs(c).labels);
            if strcmpi(p.Results.method,'bky'),
                %print uncorrected p-values for bky method since I don't
                %think there's a way to compute them
                for t=1:length(use_tpts),
                    fprintf(fid,' %f',p_values(chan_ct,t));
                end
                fprintf(fid,' p-value(uncorrected)');
            else
                for t=1:length(use_tpts),
                    fprintf(fid,' %f',adj_p(chan_ct,t));
                end
                fprintf(fid,' p-value(FDR_corrected)');
            end
            
            if mean_wind,
                %write a couple spaces and then write t-scores if mean amp
                %in time windows used
                fprintf(fid,' ');
                for t=1:n_wind
                    fprintf(fid,' %f',data_t(chan_ct,t));
                end
                fprintf(fid,' t-score \n');
            else
                fprintf(fid,'\n');
            end
        end
        
        % Write t-test details
        fprintf(fid,'Experiment: %s\n',GND.exp_desc);
        fprintf(fid,'Test_of_null_hypothesis_that_Bin_%d_equals: %d\n',bin,p.Results.null_mean);
        fprintf(fid,'#_of_time_windows: %d\n',n_wind);
        fprintf(fid,'FDR_method: %s\n',p.Results.method);
        fprintf(fid,'Tail_of_test: ');
        if ~p.Results.tail,
            fprintf(fid,'Two_tailed\n');
            fprintf(fid,'Critical_t_scores: %f %f\n',crit_t(1),crit_t(2));
        elseif p.Results.tail>0
            fprintf(fid,'Upper_tailed\n');
            fprintf(fid,'Critical_t_score: %f\n',crit_t(1));
        else
            fprintf(fid,'Lower_tailed\n');
            fprintf(fid,'Critical_t_score: %f\n',crit_t(1));
        end
        fprintf(fid,'Critical_unadjusted_p-value: %f\n',crit_p);
        fprintf(fid,'Total_#_of_significant_differences: %d\n',sum(sum(h_rej)));
        fprintf(fid,'Estimated_upper_bound_on_expected_#_of_false_discoveries: %.1f\n',sum(sum(h_rej))*p.Results.q);
        fprintf(fid,'Degrees_of_freedom: %d\n',length(use_subs)-1);
        fprintf(fid,'q_level: %f\n',p.Results.q);
        
        % # of participants and filenames
        fprintf(fid,'#_of_participants: %d\n',length(use_subs));
        fprintf(fid,'Participant_names: \n');
        for s=1:length(use_subs),
            fprintf(fid,'%s\n',GND.indiv_subnames{use_subs(s)});
        end
    end
    fclose(fid);
end

if ~strcmpi(p.Results.save_GND,'no'),
    GND=save_erps(GND);
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

