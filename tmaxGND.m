% tmaxGND() - Tests the null hypothesis that the grand average voltage
%             of a bin is mu or that the grand average within-subject 
%             difference between two bins is mu using a permutation test
%             based on the t_max statistic (e.g., Hemmelmann et al.,
%             2004).  Note, mu is assumed to be 0 by default. This function
%             requires individual subject ERPs to be stored in a  
%             "GND" structure and outputs the test results in a number of 
%             graphical and text formats.  For analogous between-subject 
%             comparisons use the function tmaxGRP.m.
%             
%             
% Usage:
%  >> [GND, prm_pval, data_t, crit_t]=tmaxGND(GND_or_fname,bin,varargin)
%
% Required Inputs:
%   GND_or_fname - A GND structure variable or the filename of a 
%                  GND structure that has been saved to disk.  To 
%                  create a GND variable from Kutaslab ERP files (e.g.,
%                  *.mas files) use avgs2GND.m.  To do the same from 
%                  EEGLAB *.set files use sets2GND.m.  See Mass
%                  Univariate ERP Toolbox documentation for detailed
%                  information about the format of a GND variable. If you
%                  specify a filename be sure to include the file's path,
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
%   alpha         - A number between 0 and 1 specifying the family-wise
%                   alpha level of the test. {default: 0.05}
%   n_perm        - The number of permutations to use in the test.  As this
%                   value increases, the test takes longer to compute and 
%                   the results become more reliable.  Manly (1997) suggests 
%                   using at least 1000 permutations for an alpha level of 
%                   0.05 and at least 5000 permutations for an alpha level 
%                   of 0.01. {default: 2500}
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
%   mean_wind     - ['yes' or 'no'] If 'yes', the permutation test will be
%                   performed on the mean amplitude within each time window 
%                   specified by time_wind.  This sacrifices temporal 
%                   resolution to increase test power by reducing the number
%                   of comparisons.  If 'no', every single time point within
%                   time_wind's time windows will be tested individually.
%                   {default: 'no'}
%   null_mean     - [number] The mean of the null hypothesis (i.e., mu) in 
%                   units of microvolts. {default: 0}
%   exclude_chans - A cell array of channel labels to exclude from the
%                   permutation test (e.g., {'A2','lle','rhe'}).  This option 
%                   sacrifices spatial resolution to increase test power by 
%                   reducing the number of comparisons. Use headinfo.m to see
%                   the channel labels stored in the GND variable. You cannot
%                   use both this option and 'include_chans' (below).{default: 
%                   not used, all channels included in test}
%   include_chans - A cell array of channel labels to use in the permutation
%                   test (e.g., {'A2','lle','rhe'}).  All other channels will
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
%                   visualizing the results of the permutation test using the 
%                   function gui_erp.m. The GUI vizualizes the grand average 
%                   ERPs in each bin via various stats (uV, t-scores), shows 
%                   topographies at individual time points, and illustrates 
%                   which electrodes significantly differ from the null 
%                   hypothesis.  This option does not work if mean_wind 
%                   option is set to 'yes.' This GUI can be reproduced using
%                   the function gui_erp.m. {default: 'yes'}
%   plot_raster   - ['yes' or 'no'] If 'yes', a two-dimensional (time x channel)
%                   binary "raster" diagram is created to illustrate the
%                   results of the permutation test.  Significant negative and
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
%                   and the details of the test (e.g., number of permutations, 
%                   family-wise alpha level, etc...). If mean_wind option is
%                   set to 'yes,' t-scores of each comparison are also 
%                   included since you cannot derive them from the t-scores
%                   at each time point/electrode in a simple way. When 
%                   importing this file into a spreadsheet be sure NOT to count
%                   consecutive spaces as multiple delimiters. {default:
%                   none}
%   save_GND      - ['yes' or 'no'] If 'yes', the GND variable will be
%                   saved to disk after the permutation test is completed 
%                   and added to it. User will first be prompted to verify 
%                   file name and path. {default: 'yes'}
%   reproduce_test- [integer] The number of the permutation test stored in
%                   the GND variable to reproduce.  For example, if 
%                   'reproduce_test' equals 2, the second t-test 
%                   stored in the GND variable (i.e., GND.t_tests(2)) will 
%                   be reproduced.  Reproduction is accomplished by setting
%                   the random number generator used in the permutation test 
%                   to the same initial state it was in when the permutation 
%                   test was first applied.
%
% Outputs:
%   GND           - GND structure variable.  This is the same as
%                   the input GND variable with one addition: the 
%                   field GND.t_tests will contain the results of the 
%                   permutation test and the test parameters. 
%   prm_pval      - A two-dimensional matrix (channel x time) of the
%                   p-values of each comparison.  These p-values are
%                   corrected for multiple comparisons by the permutation
%                   test.
%   data_t        - A two-dimensional matrix (channel x time) of the
%                   t-scores of each comparison.
%   crit_t        - The critical t-score(s) for the test.  Any t-scores
%                   that are more extreme than the critical t-score(s) 
%                   significantly deviate from 0.
%
%   Note also that a great deal of information about the test is displayed 
%   in the MATLAB command window.  You can easiy record of all this
%   information to a text file using the MATLAB command "diary."
%
% Global Variables:
%   VERBLEVEL - Mass Univariate ERP Toolbox level of verbosity (i.e., tells 
%               functions how much to report about what they're doing during
%               runtime) set by the optional function argument 'verblevel'
%
% Notes:
% -To add a difference wave to a GND variable, use the function "bin_dif".
%
% -Unlike a parametric test (e.g., an ANOVA), a discrete set of p-values
% are possible (at most the number of possible permutations).  Since the
% number of possible permutations grows rapdily with the number of
% participants, this is only issue for small sample sizes (e.g., 6
% participants).  When you have such a small sample size, the
% limited number of p-values may make the test less conservative (e.g., 
% you might be forced to use an alpha level of .0286 since it is the biggest
% possible alpha level less than .05).
%
%
% References:
%   Manly, B.F.J. (1997) Randomization, bootstrap, and Monte Carlo methods in
% biology. 2nd ed. Chapman and Hall, London.
%
%   Hemmelmann, et al. (2004) Multivariate tests for the evaluation of
% high-dimensional EEG data. Journal of Neuroscience Methods.
%
%
% Author:
% David Groppe
% Kutaslab, 3/2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 3/9/2010-'reproduce_test' option added -DG
%
% 3/14/2010-Added (optional) prompt to save GND variable at end of
% function. Added option to call sig_topo.m.  t-scores output to
% output_file along with p-values if testing mean amplitude within time
% window(s).
%
% 3/27/2010-'exclude_chans' field of GND.perm_test removed (now there's
% just a field for "include chans"
%
% 4/1/2010-'null_mean' option added.
%
% 5/5/2010-Revised to be compatible with FDR correction code 
%
% 4/13/2011-Now reports min AND max significant p-values to command line 
%



function [GND, prm_pval, data_t, crit_t]=tmaxGND(GND_or_fname,bin,varargin)

global VERBLEVEL;

p=inputParser;
p.addRequired('GND_or_fname',@(x) ischar(x) || isstruct(x));
p.addRequired('bin',@(x) isnumeric(x) && (length(x)==1) && (x>0));
p.addParamValue('tail',0,@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('alpha',0.05,@(x) isnumeric(x) && (x>0) && (x<1));
p.addParamValue('time_wind',[],@(x) isnumeric(x) && (size(x,2)==2));
p.addParamValue('time_block_dur',[],@isnumeric);
p.addParamValue('mean_wind','no',@(x) ischar(x) && (strcmpi(x,'yes') || strcmpi(x,'no')));
p.addParamValue('n_perm',2500,@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('null_mean',0,@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('verblevel',[],@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('exclude_chans',[],@(x) ischar(x) || iscell(x));
p.addParamValue('include_chans',[],@(x) ischar(x) || iscell(x));
p.addParamValue('plot_gui','yes',@(x) ischar(x) && (strcmpi(x,'yes') || strcmpi(x,'no')));
p.addParamValue('plot_raster','yes',@(x) ischar(x) && (strcmpi(x,'yes') || strcmpi(x,'no')));
p.addParamValue('plot_mn_topo',[],@(x) ischar(x) && (strcmpi(x,'yes') || strcmpi(x,'no')));
p.addParamValue('output_file',[],@ischar);
p.addParamValue('reproduce_test',[],@(x) isnumeric(x) && (length(x)==1));
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


%% Report tail of test & alpha levels
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


%% Optionally reset random number stream to reproduce a previous test
if isempty(p.Results.reproduce_test),
    seed_state=[];
else
    if p.Results.reproduce_test>length(GND.t_tests),
        error('Value of argument ''reproduce_test'' is too high.  You only have %d t-tests stored with this GND variable.',length(GND.t_tests));
    else
        if isnan(GND.t_tests(p.Results.reproduce_test).n_perm)
            error('t-test set %d is NOT a permutation test. You don''t need to seed the random number generator to reproduce it.', ...
                p.Results.reproduce_test);
        else
            seed_state=GND.t_tests(p.Results.reproduce_test).seed_state;
        end
    end
end

%Compute the permutation test
[prm_pval, data_t, crit_t, seed_state, est_alpha]=mxt_perm1(erps-p.Results.null_mean,p.Results.n_perm,p.Results.alpha,p.Results.tail,VERBLEVEL,seed_state);

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
VerbReport(sprintf('Bonferroni test-wise alpha would be %f.',p.Results.alpha/(size(prm_pval,1)* ...
    size(prm_pval,2))),1,VERBLEVEL);
sig_tpts=find(sum(prm_pval<p.Results.alpha));
if isempty(sig_tpts),
    fprintf('ERPs are NOT significantly different from zero (alpha=%f) at any time point/window analyzed.\n', ...
        p.Results.alpha);
    fprintf('All p-values>=%f\n',min(min(prm_pval)));
else
    fprintf('Significant differences from zero (in order of earliest to latest):\n');
    max_sig_p=0;
    min_sig_p=2;
    for t=sig_tpts,
        if mean_wind
            %time windows instead of time points
            fprintf('%d to %d ms, electrode(s): ',GND.time_pts(use_tpts{t}(1)), ...
                GND.time_pts(use_tpts{t}(end)));
        else
            fprintf('%d ms, electrode(s): ',GND.time_pts(use_tpts(t)));
        end
        sig_elec=find(prm_pval(:,t)<p.Results.alpha);
        ct=0;
        for c=sig_elec',
            ct=ct+1;
            if prm_pval(c,t)>max_sig_p,
                max_sig_p=prm_pval(c,t);
            end
            if prm_pval(c,t)<min_sig_p,
                min_sig_p=prm_pval(c,t);
            end
            
            if ct==length(sig_elec),
                fprintf('%s.\n',GND.chanlocs(use_chans(c)).labels);
            else
                fprintf('%s, ',GND.chanlocs(use_chans(c)).labels);
            end
        end
    end
    fprintf('All significant corrected p-values are between %f and %f\n',max_sig_p,min_sig_p);
end


%Add permutation results to GND struct
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
GND.t_tests(neo_test).mult_comp_method='tmax perm test';
GND.t_tests(neo_test).n_perm=p.Results.n_perm;
GND.t_tests(neo_test).desired_alphaORq=p.Results.alpha;
GND.t_tests(neo_test).estimated_alpha=est_alpha;
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
GND.t_tests(neo_test).adj_pval=prm_pval;
GND.t_tests(neo_test).fdr_rej=NaN;
GND.t_tests(neo_test).seed_state=seed_state;
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
    %plot_pvals(GND,p.Results.alpha,use_chans,use_tpts,prm_pval,bin); %improve plot_pvals and make an option ??
    if strcmpi(p.Results.plot_gui,'yes'),
        gui_erp('initialize','GNDorGRP',GND,'t_test',neo_test,'stat','t', ...
            'verblevel',1);
    end
end

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
            for t=1:length(use_tpts),
                fprintf(fid,' %f',prm_pval(chan_ct,t));
            end
            fprintf(fid,' p-value');
            
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
        
        % Write permutation test details
        fprintf(fid,'Experiment: %s\n',GND.exp_desc);
        fprintf(fid,'Test_of_null_hypothesis_that_Bin_%d_equals: %f\n',bin,p.Results.null_mean);
        fprintf(fid,'#_of_time_windows: %d\n',n_wind);
        fprintf(fid,'#_of_permutations: %d\n',p.Results.n_perm);
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
        fprintf(fid,'Degrees_of_freedom: %d\n',length(use_subs)-1);
        fprintf(fid,'Alpha_level: %f\n',p.Results.alpha);
        
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
    GND=save_matmk(GND);
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

