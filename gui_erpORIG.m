% gui_erp() - Open a GUI for visualzing ERP butterfly plots and
%             topographies.  ERP t-scores, standard error and global field 
%             power (Lehmann & Skrandies, 1980) can also be visualized. Hold 
%             mouse cursor over GUI controls for an explanation of what they 
%             do.  Click on waveform axes to visualize the scalp 
%             topography at that point in time. Click on electrodes in scalp 
%             topography to see electrode name.  Note, you can open more 
%             than one "gui_erp" at a time in different windows, but you 
%             cannot animate more than one gui_erp window at a time.
%
% Usage:
%  >> gui_erp(cmnd_str,varargin);
%   or
%  >> gui_erp(GND_or_GRP,varargin);
%
% Required Input:
%   cmnd_str   - One of the following strings:
%                   1.  'initialize'(intialize the GUI)
%                   2.  'time jump'(change the time point whose topography
%                        is visualized)
%                   3.  'back one'(go back one time point)
%                   4.  'forward one'(go forward one time point)
%                   5.  'back'(animate topography going backwards in time)
%                   6.  'forward'(animate topography going forwards in time)
%                   7.  'change stat'(change the ERP statistic that is being
%                        visualized [t-score or standard error])
%                   8.  'new time limits'(change the x-axis range on the
%                        waveform x time axes)
%                   9.  'new statistic limitsA'(change the y-axis range on the
%                        ERP x time axis--i.e., Axis A)
%                   10. 'new statistic limitsB' (change the y-axis range on the
%                        t-score/stderr x time axis--i.e., Axis B)
%                   11. 'change bin' (visualize a different bin)
%                   12. 'update dashed lines' (redraws dashed lines
%                        representing t-test time window and critical
%                        t-scores)
%                   13. 'redraw topo' (redraws ERP/t-score topography)
% 
%
% Optional Inputs:
%   fig_id         - [integer] ID number of the figure window in which GUI is
%                    displayed (or will be displayed).
%   bin            - [integer] The ID number of the bin whose ERPs will be
%                    visualized. Note, Kutaslab Bin 0 is ignored and is 
%                    assumed to contain cal pulses. Thus bin indexing starts 
%                    at 1.  Use headinfo.m to see the set of bins stored in 
%                    a GND or GRP variable
%   t_test         - [integer] The ID number of the set of t-tests whose
%                    results will be visualized. If specified, any optional 
%                    input arguments inconsistent with this test's parameters
%                    (e.g., 'bin') will be ignored.  Note, the results of
%                    t-tests performed on voltages averaged across 
%                    windows is not possible with this GUI; use sig_topo.m 
%                    or sig_raster.m instead. Use headinfo.m to see the 
%                    sets of t-test results stored in a GND or GRP variable      
%   show_wind      - [vector] Two numbers (in milliseconds) indicating
%                    beginning and end time points to visualize (e.g., [20 800]);
%   GNDorGRP       - MATLABmk GND structure variable containing the
%                    ERP/t-score data and t-test results to visualize.
%   critical_t     - [vector] One or two numbers indicating critical
%                    t-score(s). Time points/channels with t-scores that 
%                    exceed the critical t-score(s) are significantly 
%                    different from the mean of the null hypothesis.
%   alpha_or_q     - [number] Number between 0 and 1 indicating the family-
%                    wise alpha or FDR q level of the critical t-scores.
%   test_wind      - [vector] Pairs of numbers (in milliseconds) indicating
%                    beginning and end time points of the t-test
%                    time window(s) (e.g., [300 500]).  Multiple time windows
%                    can be specified by using semicolons to separate pairs
%                    of time points (e.g., [160 180; 300 500]).
%   ydir           - [-1 or 1] If -1, negative voltage is plotted up.  If 1 
%                    positive voltage is plotted up. {default: -1}  
%   stat           - ['t', 'gfp', or 'stder'] Statistic shown in lower GUI axis.
%                    If 't', ERPs will be shown in units of t-scores. If 'GFP',
%                    ERPs will be visualized using global field power. If 
%                    'stder', the standard error of the ERPs will be showin 
%                    in units of microvolts.
%   exclude_chans  - Cell array of channel labels to exclude from
%                    visualization (e.g., {'A2','lle','rhe'}).  If only one
%                    channel, a single string (e.g., 'A2') is acceptable.
%                    You cannot use this option AND the 't_test' option.
%   include_chans  - Cell array of channel labels to include in the
%                    visualization (e.g., {'MiPf','MiCe','MiPa','MiOc'}). If 
%                    only one channel, a single string (e.g., 'MiCe') is 
%                    acceptable.  All other channels will be ignored. You
%                    cannot use this option AND the 't_test' option.
%
% 
% Outputs:
%   None
%
%
% Author:
% David Groppe
% Kutaslab, 2/2010
%
%
% Notes:
% -If you try to use this function to visualize only channels that
% have no scalp coordinates (e.g., the difference between left and right
% hemisphere homologues), the function will automatically use sig_wave.m or
% plot_wave.m instead.
%
% -If you try to use this function to visualize only two channels that
% have scalp coordinates (e.g., MiPf & MiCe), the function will 
% automatically plot the channel locations instead of the scalp topographies.
% This is because you can't interpolate a topography with only two
% channels.  sig_wave.m is probably a better method than gui_pow.m for 
% visualizing test results at only two channels since when you click on the 
% waveforms the name of the channel corresponding to the waveform appears 
% (i.e., it's easier to determine which waveform corresponds to which channel).
% plot_wave.m is also better than gui_erp.m for one or two channels of
% data.
%
%
% References:
%    Lehmann D. & Skrandies, W. (1980) Reference-free identification of
% components of checkerboard-evoked multichannel potential fields.
% Electroencephalography and Clinical Neurophysiology. 48:609-621.



%%%%%%%%%%%%%%%%%% NOTES %%%%%%%%%%%%%%%%%%%%%
% -Much of the mechanics of this GUI is based on information that is stored
% in the figure's 'userdata' field.  The fields of the 'userdata are as
% follows:
% >>dat=get(gcf,'userdata')
% dat = 
% 
%            fig_id: 1<-The Figure ID # of the GUI
%        psbl_tests: [1 2 3 4]<-The sets of t-tests stored in the GND or
%                    GRP variable that can be visualized (mean time window
%                    tests or tests based on channels not loaded cannot be
%                    visualized)
%           t_tests: [1x4 struct]<-GND/GRP.t_test info for the possible tests 
%     mltplcmp_crct: 'fdr'<-General method that was used to correct for multiple
%                    comparisons ('fdr' or 'perm')
%             alpha: 0.0532<-q or estimated alpha level of the currently
%                    visualized test
%            n_wind: 1<-# of time windows of the currently visualized test
%        critical_t: [-9.4493 9.4493]<-critical t-scores of the currently
%                    visualized test
%         null_mean: 0<-mean of the null hypothesis of the currently
%                    visualized test
%               erp: [31x256x42 double]<-GND/GRP.grands
%          t_scores: [31x256x42 double]<-GND/GRP.grands_t
%             stder: [31x256x42 double]<-GND/GRP.grands_stder
%         showing_t: [26x256 double]<- t-scores of the data in dat.showing.
%                    These will differ from the values in dat.t_scores if 
%                    the mean of the null hypothesis of the test being 
%                    visualized is not 0
%          showingB: [26x256 double]<-waveforms currently shown in lower 
%                    butterfly plot in GUI
%          showingA: [26x256 double]<-waveforms currently shown in upper 
%                    butterfly plot in GUI
%           bindesc: {1x42 cell}<-bin descriptors for dat.erp
%             times: [1x256 double]<-time points (in ms) for dat.erp
%         plt_times: [-100 920]<-start and stop end points for waveforms
%                    currently shown
%          chanlocs: [1x31 struct]<-GND/GRP.chanlocs
%      loaded_chans: [3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
%                    23 24 25 26 27 28]<-biggest possible set of channels
%                    that can be visualized
%     showing_chans: [3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 
%                    23 24 25 26 27 28]<-set of channels currently
%                    visualized
%           h_timeA: 1.0099<-handle of upper butterfly plot axis
%          start_pt: 1<-first waveform time point to visualize
%            end_pt: 256<-last waveform time point to visualize
%            absmxA: 17.1658<-maximum value of absolute value of waveforms
%                    being visualized for upper butterfly plot
%         h_t0lineA: 3.0099<-handle for line marking time=0 in upper
%                    butterfly plot
%        h_showingA: [31x1 double]<-[26x1 double]<-handle for shown
%                    waveforms in upper butterfly plot           
%           h_lineA: 35.0099<-handle for vertical line indicating which
%                    topography is being visualized in upper butterfly plot
%          h_wind1A: 36.0099<-handle for start points of test time
%                    window(s) in upper butterfly plot  
%          h_wind2A: 37.0099<-handle for stop points of test time
%                    window(s) in upper butterfly plot  
%      h_time_ylabA: 38.0099<-handle for y-axis label on the waveform axis
%                    for upper butterfly plot  
%      h_time_title: 39.0099<-handle for title on the waveform axis
%           h_timeB: 41.0099<-handle of lower butterfly plot axis
%            absmxB: 14.5909<-maximum value of absolute value of waveforms
%                    being visualized for lower butterfly plot
%         h_t0lineB: 43.0099<-handle for line marking time=0 in lower
%                    butterfly plot
%        h_showingB: [31x1 double]<-handle for shown waveforms in lower
%                    butterfly plot 
%           h_lineB: 75.0099<-handle for vertical line indicating which
%                    topography is being visualized in lower butterfly plot
%          h_wind1B: 76.0099<-handle for start points of test time
%                    window(s) in lower butterfly plot  
%          h_wind2B: 77.0099<-handle for stop points of test time
%                    window(s) in lower butterfly plot  
%           h_crit1: 78.0099<-handle for red horizontal dashed line indicating
%                    critical t-score
%           h_crit2: 80.0099<-handle for 2nd red horizontal dashed line indicating
%                    critical t-score (used if two tailed test [i.e., two
%                    critical t-scores])
%           h_alph1: 79.0099<-handle for text box indicating alpha level.
%           h_alph2: 81.0099<-handle for 2nd text box indicating alpha level.
%                    (used if two tailed test [i.e., two critical
%                    t-scores])
%      h_time_ylabB: 83.0099<-handle for y-axis label on the waveform axis
%                    for lower butterfly plot 
%           h_topoA: 85.0099<-handle of upper topography axis
%           h_cbarA: 167.0099<-handle of topography colorbar legend for 
%                    upper axis
%       h_topo_time: 171.0099<-handle of text box indicating the time point
%                    whose topography is currently being visualized
%           h_topoB: 173.0099<-handle of lower topography axis
%           h_cbarB: 254.0099<-handle of topography colorbar legend for 
%                    lower axis
%           h_back1: 257.0099<-handle of button that moves the vertical
%                    waveform line (i.e., the point whose topography
%                    is currently being visualized) back one time point
%        h_forward1: 258.0099<-handle of button that moves the vertical
%                    waveform line (i.e., the point whose topography
%                    is currently being visualized) forward one time
%                    point
%            h_back: 259.0099<-handle of button that animates the vertical
%                    waveform line (i.e., the point whose topography
%                    is currently being visualized) back in time
%         h_forward: 260.0099<-handle of button that animates the vertical
%                    waveform line (i.e., the point whose topography
%                    is currently being visualized) forward in time
%            h_stop: 261.0099<-handle of button that stops animation of the
%                    vertical waveform line (i.e., the point whose 
%                    topography is currently being visualized)
%             h_bin: 264.0099<-handle of the scroll menu indicating which 
%                    bin is being visualized 
%           h_ptest: 266.0099<-handle of the scroll menu indicating which
%                    set of t-tests are being visualized 
%        h_testwind: 268.0099<-handle of the text box indicating the
%                    boundaries of the t-test time windows
%         h_critval: 270.0099<-handle of the text box indicating the
%                    critical t-score values of the set of t-tests
%            h_stat: 272.0099<-handle of the scroll menu indicating which
%                    statistic to plot (t-score of ERPs, standard
%                    error of ERPs)
%       h_timerange: 274.0099<-handle of the text box indicating the
%                    minimum and maximum time values on the waveform axis      
%      h_statrangeA: 276.0099<-handle of the text box indicating the
%                    minimum and maximum statistics values on the
%                    upper butterfly plot (i.e., Axis A: ERPs)
%      h_statrangeB: 278.0099<-handle of the text box indicating the
%                    minimum and maximum statistics values on the
%                    lower butterfly plot (i.e., Axis A: t-scores, 
%                    standard error)
%          help_msg: [1x447 char]<-message displayed when help button is
%                    pressed
%         interrupt: 1<- field for topography animation


%%%%%%%%% POSSIBLE FUTURE DEVELOPMENT %%%%%%%%%
%
%-Make it possible to animate more than one GUI at a time?
%
%-Make it possible to change topo color limits?  I don't know how do-able
%this is since values that exceed topo color limits are impossible to
%distinguish from values that are at the topo color limits.  Currently topo
%is automatically scaled to +/- the absolute maxima of the waveforms
%currently shown (parts of the waveform before/after the "showing time
%window" are ignored).
%
%-Make it possible to right-click on axis to get them to pop-up in a new
%figure?
%
%-Make color scheme same as EEGLAB GUIs (see variabl frm_col)?
%
%-Make time x erp axis title an optional input argument?
%
%-Make topo scale an option
%
%-Set topo color scale to green or empty when only one channel?, or delete
% it?
%
%-Add Cohen's d as a possible statistic?
%
%-Possibly fix use of axes command (see yellow warnings) to speed up
% animations?
%
%-Make frm_col=[1 1 1]*.702 the same light blue color as EEGLAB later
% (currently it's gray)?

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 4/9/2010-GUI now can visualize standard error, has a permutation test menu,
% can deal with tests on multiple time windows and tests with a non-zero 
% null mean.  Code was tidied up a bit as well.
%
% 4/12/2010-GUI now simultaneously visualizes ERPs and t-score of ERPs or
% standard error of ERPs
%
% 5/5/2010-Revised to be compatible with FDR correction code. Revised to 
% accept FDR control of t-tests. 
%
% 5/12/2010-Revised so that topographies update at same time (instead of
% staggered)
%
% 10/4/2010-Add global field power as possible statistic with which to
% visualize ERPs in the lower waveform plot.  GFP is the only waveform
% option usable if there's only one participant in the bin being visualized
%
% 10/25/2010-Revised to be able to deal with GRP variables that don't
% contain data from sufficient subjects for calculating t-scores
%
% 11/1/2010-Now checks to make sure GND or GRP variable contains at least
% one ERP and returns an error if not
%
% 12/10/2010-'show_wind' option was ineffectual. Now it works.
%
% 12/17/2010-topoplotMK can't plot topographies if there are only two
% channels of data.  Now, if there are only two channels with scalp
% coordinates, only the electrode locations are shown. Code checks to make 
% sure requested bin exists 
% 
% 3/9/2011-Code now checks to make sure there are data in a bin before
% trying to visualize them.
% 
% 6/3/2011-Fixed minor bug with viewing GRP variables and made compatible
% with cluster-based tests
% 
% 3/26/2013-Function would throw an error when trying to visualize a test
% with FDR control that didn't produce any significant results.  Should be
% fixed now.

function gui_erp(cmnd_str,varargin)

p=inputParser;
p.addRequired('cmnd_str',@(x) ischar(x) || isstruct(x));
p.addParamValue('fig_id',[],@(x) isempty(x) || (isnumeric(x) && (length(x)==1)));
p.addParamValue('bin',1,@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('show_wind',[],@(x) isempty(x) || (isnumeric(x) && (length(x)==2)));
p.addParamValue('GNDorGRP',[],@isstruct);
p.addParamValue('critical_t',[],@(x) isempty(x) || (isnumeric(x) && (length(x)<=2)));
p.addParamValue('test_wind',[],@(x) isempty(x) || (isnumeric(x) && (size(x,2)==2)));
p.addParamValue('ydir',-1,@(x) isnumeric(x) && (length(x)==1));
p.addParamValue('stat','t',@(x) ischar(x) || strcmpi(x,'t') || strcmpi(x,'stder'));
p.addParamValue('exclude_chans',[],@(x) ischar(x) || iscell(x) || isempty(x));
p.addParamValue('include_chans',[],@(x) ischar(x) || iscell(x) || isempty(x));
p.addParamValue('alpha_or_q',[],@(x) isempty(x) || isnumeric(x));
p.addParamValue('t_test',[],@(x) isempty(x) || (isnumeric(x) && (length(x)==1) && (x>0)));
p.addParamValue('verblevel',2,@(x) isnumeric(x) && (length(x)==1));
p.parse(cmnd_str,varargin{:});


if isstruct(cmnd_str)
    % Assume cmnd_str is really a GND or GRP variable and that the desired 
    % command string is 'initialize'.  Call gui_erp with same
    % optional arguments
    gui_erp('initialize','GNDorGRP',cmnd_str,'fig_id',p.Results.fig_id, ...
        'bin',p.Results.bin, ...
        'critical_t',p.Results.critical_t, ...
        'test_wind',p.Results.test_wind, ...
        'show_wind',p.Results.show_wind, ...
        'ydir',p.Results.ydir, ...
        'stat',p.Results.stat, ...
        'exclude_chans',p.Results.exclude_chans, ...
        'include_chans',p.Results.include_chans, ...
        'alpha_or_q',p.Results.alpha_or_q, ...
        't_test',p.Results.t_test, ...
        'verblevel',p.Results.verblevel);
    return;
end

% Ensure passed command is legal
psbl_cmnds{1}='initialize';
psbl_cmnds{2}='time jump';
psbl_cmnds{3}='back one';
psbl_cmnds{4}='forward one';
psbl_cmnds{5}='back';
psbl_cmnds{6}='forward';
psbl_cmnds{7}='change stat';
psbl_cmnds{8}='new time limits';
psbl_cmnds{9}='new statistic limitsA';
psbl_cmnds{10}='change bin';
psbl_cmnds{11}='new statistic limitsB';
psbl_cmnds{12}='update dashed lines';
psbl_cmnds{13}='update critical t';
psbl_cmnds{14}='redraw topo';
psbl_cmnds{15}='change test';
if ~ismember(cmnd_str,psbl_cmnds),
    error('''%s'' is not a valid value of cmnd_str for gui_erp.m',cmnd_str); 
end


if ~strcmp(cmnd_str,'initialize')
    if isempty(p.Results.fig_id),
        fig_id=gcbf;
    else
        fig_id=p.Results.fig_id;
    end
    if ~strcmp(get(fig_id,'tag'),'gui_erp')
        % If the current figure does not have the right
        % tag, find the one that does.
        h_figs = get(0,'children');
        fig_id = findobj(h_figs,'flat',...
            'tag','gui_erp');
        if isempty(fig_id),
            % If gui_erp does not exist
            % initialize it. Then run the command string
            % that was originally requested.
            gui_erp('initialize');
            gui_erp(cmnd_str);
            return;
        end
    end
    
    % At this point we know that h_fig is the handle
    % to a figure containing the GUI of interest to
    % this function (it's possible that more than one of these
    % GUI figures is open).  Therefore we can use this figure
    % handle to cut down on the number of objects
    % that need to be searched for tag names as follows:
    dat=get(gcbf,'userdata');
end

% Manage VERBLEVEL
if isempty(p.Results.verblevel),
    VERBLEVEL=2; %not global, just local
else
    VERBLEVEL=p.Results.verblevel;
end

% INITIALIZE THE GUI SECTION.
if strcmp(cmnd_str,'initialize')
    if isempty(p.Results.GNDorGRP.grands)
        error('There are no bins (i.e., ERPs) in this GND or GRP variable.  You need to add bins before you can visualize them with the ERP GUI.');
    end
    
    % Creates a new GUI (even if one already exists)
    if isempty(p.Results.fig_id),
        dat.fig_id=figure;
    else
        dat.fig_id=p.Results.fig_id;
        figure(dat.fig_id);
        clf
    end
    set(dat.fig_id,'name',['ERP GUI ' p.Results.GNDorGRP.exp_desc],'tag','gui_erp', ...
        'MenuBar','none','position',[139 42 690 700]);
    
    %% Manage t-test results
    n_t_tests=length(p.Results.GNDorGRP.t_tests);
    dat.psbl_tests=[];
    temp_bins=zeros(1,n_t_tests);
    for d=1:n_t_tests,
        if ~strcmpi(p.Results.GNDorGRP.t_tests(d).mean_wind,'yes'),
            dat.psbl_tests=[dat.psbl_tests d];
            temp_bins(d)=p.Results.GNDorGRP.t_tests(d).bin;
        end
    end
    temp_bins=temp_bins(dat.psbl_tests);
    n_psbl_tests=length(dat.psbl_tests);
    dat.t_tests=p.Results.GNDorGRP.t_tests(dat.psbl_tests);
    crnt_ttest=n_psbl_tests+1; %set of t-tests being visualized, default is None/Manual
    
    use_ttest=0;
    if isempty(p.Results.t_test),
        if isempty(p.Results.include_chans) && isempty(p.Results.exclude_chans) ...
                && isempty(p.Results.critical_t) && isempty(p.Results.alpha_or_q) ...
                && isempty(p.Results.test_wind),
            %search for set of t-tests that have been performed on bin being
            %visualized (if potentially incompatible optional input
            %arguments have NOT been specified)
            test_ids=find(p.Results.bin==temp_bins);
            if ~isempty(test_ids),
                use_ttest=1;
                crnt_ttest=test_ids(1); %in case there's more than one test
                if VERBLEVEL>=2
                    fprintf('Plotting results of t-tests set %d.  Any input arguments (e.g., ''critical_t'') inconsistent with this test will be ignored.\n', ...
                        dat.psbl_tests(crnt_ttest));
                end
            end
        end
    else
        if p.Results.t_test>n_t_tests,
            error('Argument ''t_test'' value of %d exceeds the number of test results stored in this GND/GRP variable (i.e., %d).', ...
                p.Results.t_test,n_t_tests);
        elseif ~ismember(p.Results.t_test,dat.psbl_tests),
            error('t-test set %d in this GND/GRP variable was performed on mean amplitudes within one or more time windows.  You cannot visualize such tests with gui_erp.m.', ...
                p.Results.t_test);
        else
            use_ttest=1;
            crnt_ttest=find(dat.psbl_tests==p.Results.t_test);
            if VERBLEVEL>=2
                fprintf('Plotting results of t-test set %d.  Any input arguments (e.g., ''bin'') inconsistent with this test will be ignored.\n', ...
                    p.Results.t_test);
            end
        end
    end
    if use_ttest,
        % set all other optional arguments to be
        % consistent with t-test 
        bin=dat.t_tests(crnt_ttest).bin;
        if isnan(dat.t_tests(crnt_ttest).estimated_alpha)
            dat.alpha=dat.t_tests(crnt_ttest).desired_alphaORq;
            dat.mltplcmp_crct='fdr';
        else
            dat.alpha=dat.t_tests(crnt_ttest).estimated_alpha;
            dat.mltplcmp_crct='perm';
        end
        test_wind=dat.t_tests(crnt_ttest).time_wind;
        dat.n_wind=size(test_wind,1);
        include_chans=dat.t_tests(crnt_ttest).include_chans;
        dat.critical_t=dat.t_tests(crnt_ttest).crit_t;
        dat.null_mean=dat.t_tests(crnt_ttest).null_mean;
    else
        % no t-test, use optional inputs or defaults
        bin=p.Results.bin;
        test_wind=p.Results.test_wind;
        include_chans=p.Results.include_chans;
        dat.alpha=p.Results.alpha_or_q;
        dat.critical_t=p.Results.critical_t;
        dat.null_mean=0;
    end
    
    
    %% Make sure user hasn't asked for a bin that doesn't exist or that
    % doesn't contain any data
    if bin>length(p.Results.GNDorGRP.bin_info)
        close(gcf);
        error('You asked to visualize Bin %d, but your GND/GRP variable only contains %d bins.', ...
            bin,length(p.Results.GNDorGRP.bin_info));
    elseif (isfield(p.Results.GNDorGRP,'sub_ct') && ~p.Results.GNDorGRP.sub_ct(bin)) || ...
            (isfield(p.Results.GNDorGRP.bin_info(bin),'n_subsA') && ( ~p.Results.GNDorGRP.bin_info(bin).n_subsA || ~p.Results.GNDorGRP.bin_info(bin).n_subsB))
        close(gcf);
        error('You asked to visualize Bin %d, but your GND/GRP variable doesn''t have any data in that bin.', ...
            bin);
    end
        
    %% Figure out which channels to ignore if any
    n_chan=length(p.Results.GNDorGRP.chanlocs);
    %Make sure exclude & include options were not both used.
    if ~isempty(include_chans) && ~isempty(p.Results.exclude_chans)
        if use_ttest,
            error('You cannot specify a set of t-tests to visualize and use the ''exclude_chans'' option.');
        else
            error('You cannot use BOTH ''include_chans'' and ''exclude_chans'' options.');
        end
    end
    if ischar(p.Results.exclude_chans),
        exclude_chans{1}=p.Results.exclude_chans;
    elseif isempty(p.Results.exclude_chans)
        exclude_chans=[];
    else
        exclude_chans=p.Results.exclude_chans;
    end
    if ischar(include_chans),
        temp_var=include_chans;
        clear include_chans;
        include_chans{1}=temp_var;
        clear temp_var;
    end
    if ~isempty(exclude_chans),
        ignore_chans=zeros(1,length(exclude_chans)); %preallocate mem
        ct=0;
        for x=1:length( exclude_chans),
            found=0;
            for c=1:n_chan,
                if strcmpi(exclude_chans{x},p.Results.GNDorGRP.chanlocs(c).labels),
                    found=1;
                    ct=ct+1;
                    ignore_chans(ct)=c;
                end
            end
            if ~found,
                watchit(sprintf('I attempted to exclude %s.  However no such electrode was found in GND/GRP variable.', ...
                    exclude_chans{x}));
            end
        end
        ignore_chans=ignore_chans(1:ct);
        loaded_chans=setdiff(1:n_chan,ignore_chans);
    elseif ~isempty(include_chans),
        loaded_chans=zeros(1,length(include_chans)); %preallocate mem
        ct=0;
        for x=1:length(include_chans),
            found=0;
            for c=1:n_chan,
                if strcmpi(include_chans{x},p.Results.GNDorGRP.chanlocs(c).labels),
                    found=1;
                    ct=ct+1;
                    loaded_chans(ct)=c;
                end
            end
            if ~found,
                watchit(sprintf('I attempted to include %s.  However no such electrode was found in GND/GRP variable.', ...
                     include_chans{x}));
            end
        end
        loaded_chans=loaded_chans(1:ct);
    else
        loaded_chans=1:n_chan;
    end
    if isempty(loaded_chans),
       error('No channels selected for visualization!'); 
    end
    %Check to see if any channels are missing coordinates
    yes_coord=zeros(1,n_chan);
    for c=loaded_chans,
        if ~isempty(p.Results.GNDorGRP.chanlocs(c).theta) && ~isnan(p.Results.GNDorGRP.chanlocs(c).theta)
            yes_coord(c)=1;
        end
    end
    if (sum(yes_coord)==0)
        if use_ttest
            %t-test results should be shown
            watchit(sprintf('None of the channels you wish to visualize have scalp coordinates.\nIt is pointless to use the ERP GUI.  Using sig_wave.m instead.'));
            close(gcf);
            sig_wave(p.Results.GNDorGRP,crnt_ttest,'ydir',p.Results.ydir, ...
                'verblevel',p.Results.verblevel);
        else
            watchit('None of the channels you wish to visualize have scalp coordinates.  It is pointless to use the ERP GUI.  Using plot_wave.m instead.');
            close(gcf);
            if ~isempty(include_chans),
                plot_wave(p.Results.GNDorGRP,bin,'ydir',p.Results.ydir, ...
                    'verblevel',p.Results.verblevel,'include_chans',include_chans);
            else
                %Note, exclude_chans should be empty if all channels are to
                %be shown
                plot_wave(p.Results.GNDorGRP,bin,'ydir',p.Results.ydir, ...
                    'verblevel',p.Results.verblevel,'exclude_chans',exclude_chans);
            end
        end
        return
    elseif sum(yes_coord)<length(loaded_chans),
        watchit(sprintf('%d channels do not have scalp coordinates.  The waveforms for such channels will be visualized but they the will not be represented in the scalp topography.', ...
            length(loaded_chans)-sum(yes_coord)));
    end
        
    %remove any sets of t-tests peformed on channels that will not be
    %included in the GUI
    use_tests=[];
    for d=1:n_psbl_tests,
       if ~isempty(intersect(loaded_chans,dat.t_tests(d).used_chan_ids)),
          use_tests=[use_tests d]; 
       end
    end
    dat.psbl_tests=dat.psbl_tests(use_tests);
    dat.t_tests=dat.t_tests(use_tests);
    crnt_ttest=find(crnt_ttest==[use_tests n_psbl_tests+1]); % "n_psbl_tests+1" adds the manual/no test option
    n_psbl_tests=length(use_tests);
    
    %% attach data to figure
    dat.erp=p.Results.GNDorGRP.grands;
    dat.t_scores=p.Results.GNDorGRP.grands_t;
    dat.stder=p.Results.GNDorGRP.grands_stder;
    mn_erp_across_chans=mean(dat.erp,1);
    dat.gfp=squeeze(sqrt(mean( (dat.erp-repmat(mn_erp_across_chans,[n_chan 1 1])).^2,1)));
    %t-scores of data that will be shown
    if (crnt_ttest<=n_psbl_tests) && (dat.t_tests(crnt_ttest).null_mean),
        %The mean of the null hypothesis is non-zero
        dat.showing_t=squeeze( (p.Results.GNDorGRP.grands(loaded_chans,:,bin)- ...
            dat.t_tests(crnt_ttest).null_mean)./p.Results.GNDorGRP.grands_stder(loaded_chans,:,bin) ); 
    else
        dat.showing_t=squeeze(p.Results.GNDorGRP.grands_t(loaded_chans,:,bin)); 
    end
    %data that will be shown
    if isinf(dat.t_scores(1,1,bin)) || isnan(dat.t_scores(1,1,bin)) || strcmpi(p.Results.stat,'gfp'),
        %default to global field power if there's only one subject in the
        %bin (and standard deviation isn't defined). For within-subject
        %t-scores, stdev will bin Inf/-Inf when there's only one subject.
        %For between-subject t-scores, stdev will bin NaN when there are
        %insufficient subjects in one or both groups.
        dat.showingB=dat.gfp(:,bin)'; 
    elseif strcmpi(p.Results.stat,'t'),
        dat.showingB=dat.showing_t;
    else
        %standard error
        dat.showingB=squeeze(p.Results.GNDorGRP.grands_stder(loaded_chans,:,bin)); 
    end
    dat.showingA=squeeze(p.Results.GNDorGRP.grands(loaded_chans,:,bin)); %always show ERPs in top axis
    
    n_bin=length(p.Results.GNDorGRP.bin_info);
    dat.bindesc=cell(1,n_bin);
    for b=1:n_bin,
        dat.bindesc{b}=p.Results.GNDorGRP.bin_info(b).bindesc;
    end
    dat.times=p.Results.GNDorGRP.time_pts;
    if isempty(p.Results.show_wind),
        dat.plt_times=[p.Results.GNDorGRP.time_pts(1) p.Results.GNDorGRP.time_pts(end)];
    else
        dat.plt_times=p.Results.show_wind;
    end
    dat.chanlocs=p.Results.GNDorGRP.chanlocs;
    dat.loaded_chans=loaded_chans;
    dat.showing_chans=loaded_chans;
    critical_t=dat.critical_t;
    ydir=p.Results.ydir; 
    if length(critical_t)>2,
       critical_t=critical_t(1:2); 
    end
    dat.critical_t=critical_t;
    
    %
    %%%%%%%%%%%%%%%%%%%%%% AXIS A: Time x ERP Axes %%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    frm_col=[1 1 1]*.702;
    uipanel(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0 0.59 .72 0.42 ],...
        'shadowcolor','k', ...
        'highlightcolor',frm_col, ...
        'foregroundcolor',frm_col, ...
        'backgroundcolor',frm_col);
    
    dat.h_timeA=axes('position',[0.1 .62 .6 .33]);
    dat.start_pt=find_tpt(dat.plt_times(1),dat.times);
    dat.end_pt=find_tpt(dat.plt_times(2),dat.times);
    plotted_pts=dat.start_pt:dat.end_pt;
    plotted_times=dat.times(dat.start_pt:dat.end_pt);
    [dat.absmxA mx_tpt]=max(max(abs(dat.showingA(:,dat.start_pt:dat.end_pt))));
    stat_mx=max(max(dat.showingA(:,dat.start_pt:dat.end_pt)));
    stat_mn=min(min(dat.showingA(:,dat.start_pt:dat.end_pt)));
    stat_rng=stat_mx-stat_mn;
    stat_plt_rng=[stat_mn-stat_rng*.02 stat_mx+stat_rng*.02];
    stat_plt_rng=round(stat_plt_rng*100)/100;
    plot([p.Results.GNDorGRP.time_pts(1) p.Results.GNDorGRP.time_pts(end)],[0 0],'k'); % uV/t=0 line
    hold on;
    dat.h_t0lineA=plot([0 0],stat_plt_rng,'k'); % time=0 line
    set(dat.h_timeA,'ygrid','on');
    if ydir<0,
       set(dat.h_timeA,'ydir','reverse'); 
    end
    if size(dat.showingA,1)==1,
        %only one channel being plot
        dat.h_showingA=plot(dat.times,dat.showingA);
    else
        dat.h_showingA=plot(dat.times,dat.showingA');
    end
    axis tight;
    axis([plotted_times(1) plotted_times(end) stat_plt_rng]);
    vA=axis;
    dat.h_lineA=plot([1 1]*plotted_times(mx_tpt),vA(3:4),'k');
    set(dat.h_lineA,'linewidth',2);
    
    % Dashed vertical lines marking test window
    dat.n_wind=size(test_wind,1);
    if dat.n_wind        
        for nw=1:dat.n_wind,
            dat.h_wind1A(nw)=plot([1 1]*test_wind(nw,1),vA(3:4),'k--');
            set(dat.h_wind1A(nw),'linewidth',2);
            dat.h_wind2A(nw)=plot([1 1]*test_wind(nw,2),vA(3:4),'k--');
            set(dat.h_wind2A(nw),'linewidth',2);
        end
    else
        dat.h_wind1A=[];
        dat.h_wind2A=[];
    end
    
    %y-axis
    h=ylabel('\muV (ERP)');
    dat.h_time_ylabA=h;
    set(h,'fontsize',10,'fontunits','normalized');
    
    %Title
    new_title=['Bin ' int2str(bin) ': ' dat.bindesc{bin}];
    title_max_char=43; %prevents title from spilling over into topography axis
    if length(new_title)>title_max_char,
        new_title=new_title(1:title_max_char);
    end 
    dat.h_time_title=title(new_title);
    set(dat.h_time_title,'fontsize',10,'fontunits','normalized');
    bdf_code = [ 'tmppos = get(gca, ''currentpoint'');' ...
        'dat=get(gcbf, ''userdata'');' ...
        'new_tpt=find_tpt(tmppos(1,1),dat.times);' ...
        'set(dat.h_lineA,''XData'',[1 1]*dat.times(new_tpt));' ...
        'set(dat.h_lineB,''XData'',[1 1]*dat.times(new_tpt));' ...
        'set(dat.h_topo_time,''string'',num2str(dat.times(new_tpt)));' ...
        'set(dat.fig_id,''userdata'',dat);' ...
        'gui_erp(''redraw topo'');' ...
        'drawnow;' ...
        'clear latpoint dattmp tmppos;' ...
        ];
    set(dat.h_timeA,'ButtonDownFcn',bdf_code);
    set(dat.h_showingA,'ButtonDownFcn',bdf_code);

    
    %
    %%%%%%%%%%%%%%%%%%%%%% AXIS B: Time x t-score/stderr/GFP Axes %%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    frm_col=[1 1 1]*.702; %background color (gray) differs from EEGLAB blue background color
    uipanel(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0 0.175 .72 0.42 ],...
        'shadowcolor','k', ...
        'highlightcolor',frm_col, ...
        'foregroundcolor',frm_col, ...
        'backgroundcolor',frm_col);
   
    dat.h_timeB=axes('position',[0.1 .25 .6 .33]); 
    dat.absmxB=max(max(abs(dat.showingB(:,dat.start_pt:dat.end_pt))));
    stat_mx=max(max(dat.showingB(:,dat.start_pt:dat.end_pt)));
    stat_mn=min(min(dat.showingB(:,dat.start_pt:dat.end_pt)));
    stat_rng=stat_mx-stat_mn;
    stat_plt_rng=[stat_mn-stat_rng*.02 stat_mx+stat_rng*.02];
    stat_plt_rng=round(stat_plt_rng*100)/100;
    plot([p.Results.GNDorGRP.time_pts(1) p.Results.GNDorGRP.time_pts(end)],[0 0],'k'); % uV/t=0 line
    hold on;
    dat.h_t0lineB=plot([0 0],stat_plt_rng,'k'); % time=0 line
    set(dat.h_timeB,'ygrid','on');
    if ydir<0,
       set(dat.h_timeB,'ydir','reverse'); 
    end
    if size(dat.showingB,1)==1,
        %only one channel being plot
        dat.h_showingB=plot(dat.times,dat.showingB);
    else
        dat.h_showingB=plot(dat.times,dat.showingB');
    end
    axis tight;
    axis([plotted_times(1) plotted_times(end) stat_plt_rng]);
    vB=axis;
    dat.h_lineB=plot([1 1]*plotted_times(mx_tpt),vB(3:4),'k');
    set(dat.h_lineB,'linewidth',2);
    
    % Dashed vertical lines marking test window
    if dat.n_wind        
        for nw=1:dat.n_wind,
            dat.h_wind1B(nw)=plot([1 1]*test_wind(nw,1),vB(3:4),'k--');
            set(dat.h_wind1B(nw),'linewidth',2);
            dat.h_wind2B(nw)=plot([1 1]*test_wind(nw,2),vB(3:4),'k--');
            set(dat.h_wind2B(nw),'linewidth',2);
        end
    else
        dat.h_wind1B=[];
        dat.h_wind2B=[];
    end
    
    % Dashed lines marking critical t-score(s)
    dat.h_crit1=[];
    dat.h_crit2=[];
    dat.h_alph1=[];
    dat.h_alph2=[];
    if ~isempty(critical_t) && ~isempty(test_wind),
        if strcmpi(p.Results.stat,'t'),
            if isempty(dat.alpha),
                watchit('You did not specify an alpha level via optional input argument ''alpha''.  Alpha level will not be displayed.');
            else
                for nw=1:dat.n_wind,
                    dat.h_crit1(nw)=plot(test_wind(nw,:),[1 1]*critical_t(1),'r--');
                    set(dat.h_crit1(nw),'linewidth',3);
                    if nw==1,
                        tm_rng=p.Results.GNDorGRP.time_pts(end)-p.Results.GNDorGRP.time_pts(1);
                        if strcmpi(dat.mltplcmp_crct,'fdr')
                            dat.h_alph1=text(test_wind(nw,1)-tm_rng*.02,critical_t(1), ...
                                ['q=' num2str(rnd_orderofmag(dat.alpha))]);
                        else
                            dat.h_alph1=text(test_wind(nw,1)-tm_rng*.02,critical_t(1), ...
                                ['\alpha=' num2str(rnd_orderofmag(dat.alpha))]);
                        end
                        set(dat.h_alph1,'color','r','fontweight','normal','fontsize',12, ...
                            'horizontalalignment','right','backgroundcolor',[1 1 1], ...
                            'edgecolor',[1 1 1]*.3,'clipping','on','fontname','fixedwidth');
                    end
                end
            end
            if length(critical_t)>1,
                for nw=1:dat.n_wind,
                    dat.h_crit2(nw)=plot(test_wind(nw,:),[1 1]*critical_t(2),'r--');
                    set(dat.h_crit2(nw),'linewidth',3);
                    if nw==1,
                        if strcmpi(dat.mltplcmp_crct,'fdr')
                            dat.h_alph2=text(test_wind(nw,1)-tm_rng*.02,critical_t(2), ...
                                ['q=' num2str(rnd_orderofmag(dat.alpha))]);
                        else
                            dat.h_alph2=text(test_wind(nw,1)-tm_rng*.02,critical_t(2), ...
                                ['\alpha=' num2str(rnd_orderofmag(dat.alpha))]);
                        end
                        set(dat.h_alph2,'color','r','fontweight','normal','fontsize',12, ...
                            'horizontalalignment','right','backgroundcolor',[1 1 1], ...
                            'edgecolor',[1 1 1]*.3,'clipping','on','fontname','fixedwidth');
                    end
                end
            end
        end
    end
    %x-axis label
    h=xlabel('Time (msec)');
    set(h,'fontsize',10,'fontunits','normalized');

    %y-axis label
    if strcmpi(p.Results.stat,'gfp') || isinf(dat.t_scores(1,1,bin)) || isnan(dat.t_scores(1,1,bin)),
        h=ylabel('\muV (GFP)');
    elseif strcmpi(p.Results.stat,'t'),
        h=ylabel('t-score');
    else
        h=ylabel('\muV (StdEr)');
    end
    dat.h_time_ylabB=h;
    set(h,'fontsize',10,'fontunits','normalized');
    
    bdf_code = [ 'tmppos = get(gca, ''currentpoint'');' ...
        'dat=get(gcbf, ''userdata'');' ...
        'new_tpt=find_tpt(tmppos(1,1),dat.times);' ...
        'set(dat.h_lineA,''XData'',[1 1]*dat.times(new_tpt));' ...
        'set(dat.h_lineB,''XData'',[1 1]*dat.times(new_tpt));' ...
        'set(dat.h_topo_time,''string'',num2str(dat.times(new_tpt)));' ...
        'set(dat.fig_id,''userdata'',dat);' ...
        'gui_erp(''redraw topo'');' ...
        'drawnow;' ...
        'clear latpoint dattmp tmppos;' ...
        ];
    set(dat.h_timeB,'ButtonDownFcn',bdf_code);
    set(dat.h_showingB,'ButtonDownFcn',bdf_code);

    %
    %%%%%%%%%%%%%%%%%%%%%% AXIS C: ERP Topography %%%%%%%%%%%%%%%%%%%%%%%%
    %

    % PANEL
    uipanel(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0.719 0.59 0.281 0.42 ],...
        'shadowcolor','k', ...
        'highlightcolor',frm_col, ...
        'foregroundcolor',frm_col, ...
        'backgroundcolor',frm_col);   
    
    % TOPOGRAPHY
    dat.h_topoA=axes('position',[0.705 .68 .31 .24],'box','off'); 
  
    sig_chans=[];
    if ~isempty(test_wind) && ~isempty(critical_t),
        %if current time point is in a test window look for channels with
        %sig effects
        for nw=1:dat.n_wind,
            if (dat.times(plotted_pts(mx_tpt))>=test_wind(nw,1)) && (dat.times(plotted_pts(mx_tpt))<=test_wind(nw,2)),
                if ~strcmpi(dat.mltplcmp_crct,'fdr') && (length(critical_t)==1) && isnan(critical_t) % critical_t is always the scalar NaN for cluster based permutation tests
                    %cluster based test
                    mx_tpt_in_ms=plotted_times(mx_tpt);
                    mx_tpt_epoch_id=find(dat.times==mx_tpt_in_ms);
                    use_test_id=find(dat.psbl_tests==p.Results.t_test);
                    pval_tpt_id=find(dat.t_tests(use_test_id).used_tpt_ids==mx_tpt_epoch_id);
                    if ~isempty(pval_tpt_id)
                        sig_chans_temp=find(dat.t_tests(use_test_id).adj_pval(:,pval_tpt_id)<dat.t_tests(use_test_id).desired_alphaORq);
                        sig_chans_temp=dat.t_tests(use_test_id).used_chan_ids(sig_chans_temp); %convert sig channel indices into channel indices in the original GND/GRP variable
                        n_showing_chans=length(dat.showing_chans);
                        sig_chans=zeros(1,n_showing_chans);
                        for a=1:n_showing_chans,
                            if ismember(dat.showing_chans(a),sig_chans_temp)
                                sig_chans(a)=1;
                            end
                        end
                        sig_chans=find(sig_chans);
                    else
                        sig_chans=[];
                    end
                else
                    if length(critical_t)==2,
                        % I don't think one critical t value can be NaN and not
                        % the other, but just in case
                        if isnan(critical_t(2))
                            sig_chans=[];
                        else
                            sig_chans=find(dat.showing_t(:,plotted_pts(mx_tpt))>max(critical_t));
                        end
                        if ~isnan(critical_t(1))
                            sig_chans=[sig_chans; find(dat.showing_t(:,plotted_pts(mx_tpt))<min(critical_t))];
                        end
                    else
                        if isnan(critical_t) %FDR control con produce NaN critical_t's if not test are significant
                            sig_chans=[];
                        elseif critical_t>0,
                            sig_chans=find(dat.showing_t(:,plotted_pts(mx_tpt))>critical_t);
                        else
                            sig_chans=find(dat.showing_t(:,plotted_pts(mx_tpt))<critical_t);
                        end
                    end
                end
                break; %Visualized time point is in this window. Break out of for loop since there's no need to look at additional time windows 
            end
        end
    end
    cbar_title_fontsize=14;
    if size(dat.showingA,1)<=2,
        %two or fewer channels, we can only plot electrode locations
        topoplotMK(dat.showingA(:,plotted_pts(mx_tpt)),dat.chanlocs(dat.showing_chans), ...
            'style','blank','plain_blank',1,'emarker2',{sig_chans,'o',[1 1 1],4});
        cbar_title='Not Applicable';
        cbar_title_fontsize=10; %for some reason a fontsize of 14 cuts off last later for topoB axis
    else
        topoplotMK(dat.showingA(:,plotted_pts(mx_tpt)),dat.chanlocs(dat.showing_chans), ...
            'maplimits',[-1 1]*dat.absmxA,'emarker2',{sig_chans,'o',[1 1 1],4});
        set(findobj(gca,'type','patch'),'facecolor',[1 1 1]*.702);
        cbar_title='\muV';
    end

    
    % COLOR BAR
    dat.h_cbarA=axes('position',[0.76 .945 .2 .015]);
    cbar(dat.h_cbarA);
    absmx=round(dat.absmxA*100)/100;
    set(gca,'xticklabel',[-absmx 0 absmx]);
    h_cbar_title=title(cbar_title);
    set(h_cbar_title,'fontsize',cbar_title_fontsize);
    
    
    %%% TIME POINT TEXT BOX
    %STATIC TEXT LABEL
    uicontrol(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0.905 0.62 0.07 0.04 ],...
        'String','msec',...
        'fontsize',14, ...
        'fontunits','normalized', ...
        'Style','text');
    dat.h_topo_time=uicontrol(dat.fig_id,...
        'CallBack','gui_erp(''time jump'');',...
        'Units','normalized', ...
        'Position',[ 0.805 0.62 0.1 0.05 ], ...
        'String',num2str(plotted_times(mx_tpt)), ...
        'Style','edit', ...
        'Enable','on', ...
        'ToolTipString','Time point to visualize topographicaly.', ...
        'fontsize',14, ...
        'fontunits','normalized', ...
        'horizontalalignment','center', ...
        'BackGroundColor','w', ...
        'Tag','topo_time');
    
    
    %
    %%%%%%%%%%%%%%%%%%%%%% AXIS D: t-Score/Standard Error/GFP Topography %%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    
    % PANEL
    uipanel(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0.719 0.175 0.281 0.42 ],...
        'shadowcolor','k', ...
        'highlightcolor',frm_col, ...
        'foregroundcolor',frm_col, ...
        'backgroundcolor',frm_col);
    
    % TOPOGRAPHY
    dat.h_topoB=axes('position',[0.705 .27 .31 .24],'box','off');
    
    cbar_title_fontsize=14;
    if strcmpi(p.Results.stat,'gfp') || isinf(dat.t_scores(1,1,bin)) || isnan(dat.t_scores(1,1,bin))
        topoplotMK([],dat.chanlocs(dat.showing_chans), ...
            'style','blank','plain_blank',1);
        cbar_title='Not Applicable';
        cbar_title_fontsize=10;
    else
        if size(dat.showingB,1)<=2,
            %two or fewer channels, we can only plot electrode locations
            topoplotMK(dat.showingB(:,plotted_pts(mx_tpt)),dat.chanlocs(dat.showing_chans), ...
                'style','blank','plain_blank',1,'emarker2',{sig_chans,'o',[1 1 1],4});
            cbar_title='Not Applicable';
            cbar_title_fontsize=10;
        else
            topoplotMK(dat.showingB(:,plotted_pts(mx_tpt)),dat.chanlocs(dat.showing_chans), ...
                'maplimits',[-1 1]*dat.absmxB,'emarker2',{sig_chans,'o',[1 1 1],4});
            set(findobj(gca,'type','patch'),'facecolor',[1 1 1]*.702);
            if strcmpi(p.Results.stat,'t'),
                cbar_title='t-score';
            else
                cbar_title='\muV'; %ERPs and stder are in units of uV
            end
        end
    end
    
    % COLOR BAR
    dat.h_cbarB=axes('position',[0.76 .535 .2 .015]);
    cbar(dat.h_cbarB);
    absmx=round(dat.absmxB*100)/100;
    set(gca,'xticklabel',[-absmx 0 absmx]);
    h_cbar_title=title(cbar_title);
    set(h_cbar_title,'fontsize',cbar_title_fontsize);
    
    % GO BACK 1 TIME POINT
    dat.h_back1=uicontrol(dat.fig_id,...
        'CallBack','gui_erp(''back one'');',...
        'Units','normalized', ...
        'Position',[ 0.80 0.232 0.055 0.04 ],...
        'String','<',...
        'fontsize',14, ...
        'fontunits','normalized', ...
        'Tag','back1', ...
        'ToolTipString','Show topography at preceding time point.', ...
        'Style','pushbutton');
    
    % GO FORWARD 1 TIME POINT
    dat.h_forward1=uicontrol(dat.fig_id,...
        'CallBack','gui_erp(''forward one'');',...
        'Units','normalized', ...
        'Position',[ 0.865 0.232 0.055 0.04 ],...
        'String','>', ...
        'fontsize',14, ...
        'fontunits','normalized', ...
        'Tag','forward1', ...
        'ToolTipString','Show topography at subsequent time point.', ...
        'Style','pushbutton');
    
    % ANIMATE BACKWARDS IN TIME
    dat.h_back=uicontrol(dat.fig_id,...
        'CallBack','gui_erp(''back'');',...
        'Units','normalized', ...
        'Position',[ 0.735 0.232 0.055 0.04 ],...
        'String','<<',...
        'fontsize',14, ...
        'fontunits','normalized', ...
        'Tag','back', ...
        'interruptible','on', ...
        'ToolTipString','Animate topography going backwards in time.', ...
        'Style','pushbutton');
    
    % ANIMATE FORWARDS IN TIME
    dat.h_forward=uicontrol(dat.fig_id,...
        'CallBack','gui_erp(''forward'');',...
        'Units','normalized', ...
        'Position',[ 0.93 0.232 0.055 0.04 ],...
        'String','>>',...
        'fontsize',14, ...
        'fontunits','normalized', ...
        'Tag','back', ...
        'interruptible','on', ...
        'ToolTipString','Animate topography going forwards in time.', ...
        'Style','pushbutton');
    
    % STOP ANIMATION
    dat.h_stop=uicontrol(dat.fig_id,...
        'CallBack','dat=get(gcbf,''userdata''); dat.interrupt=1; set(dat.fig_id,''userdata'',dat);',...
        'Units','normalized', ...
        'Position',[ 0.745 0.179 0.23 0.05 ],...
        'String','Stop Animation',...
        'fontsize',14, ...
        'fontunits','normalized', ...
        'Tag','stop', ...
        'enable','off', ...
        'ToolTipString','Stop topography animation.', ...
        'Style','pushbutton');
   

    
    %
    %%%%%%%%%%%%%%%%%%%%%% AXIS E: PLOTTING OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    % PANEL
    uipanel(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0 0 1 0.178 ],...
        'shadowcolor','k', ...
        'highlightcolor',frm_col, ...
        'foregroundcolor',frm_col, ...
        'backgroundcolor',frm_col);
      
        
    %%%% WHICH BIN TO PLOT
    bdesc_str=cell(1,n_bin);
    for a=1:n_bin,
        bdesc_str{a}=sprintf('Bin %d: %s',a,dat.bindesc{a});
    end 
    % Text
    uicontrol(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0.01 0.115 0.2 0.05 ],...
        'String','Current Bin',...
        'fontsize',12, ...
        'ToolTipString','Bin to plot in Time x ERP/t-score axis.', ...
        'Style','text');
    % Create bin menu
    dat.h_bin=uicontrol(dat.fig_id,...
        'CallBack','gui_erp(''change bin'');',...
        'Units','normalized', ...
        'Position',[ 0.01 0.08 0.24 0.05 ],...
        'fontsize',12, ...
        'String',bdesc_str, ...
        'Value',bin, ...
        'Style','popup', ...
        'ToolTipString','Bin to plot in Time x ERP/t-score axis.', ...
        'tag','bin');
    
    
    %%%% WHICH T-TESTS RESULTS TO PLOT
    ttest_str=cell(1,n_psbl_tests+1);
    for a=1:n_psbl_tests,
        ttest_str{a}=sprintf('Test %d: Bin %d, Times ',dat.psbl_tests(a),dat.t_tests(a).bin);
        n_ttst_wind=size(dat.t_tests(a).time_wind,1);
        for b=1:n_ttst_wind,
            ttest_str{a}=[ttest_str{a} int2str(dat.t_tests(a).time_wind(b,1)) '-' int2str(dat.t_tests(a).time_wind(b,2)) ','];
        end
        if length(dat.t_tests(a).crit_t)==2,
            ttest_str{a}=[ttest_str{a} ' Two-Tailed,'];
        elseif dat.t_tests(a).crit_t<0,
            ttest_str{a}=[ttest_str{a} ' Lower-Tailed,'];
        else
            ttest_str{a}=[ttest_str{a} ' Upper-Tailed,'];
        end
        
        if isnan(dat.t_tests(a).estimated_alpha)
            ttest_str{a}=[ttest_str{a} sprintf(' q=%.4f',dat.t_tests(a).desired_alphaORq)];
            ttest_str{a}=[ttest_str{a} sprintf(', method=%s',dat.t_tests(a).mult_comp_method)];
        else
            ttest_str{a}=[ttest_str{a} sprintf(' alpha=%.4f',dat.t_tests(a).estimated_alpha)];
        end
        
        if dat.t_tests(a).null_mean,
            ttest_str{a}=[ttest_str{a} sprintf(', Null Mean=%.4f',dat.t_tests(a).null_mean)];
        end
    end
    ttest_str{n_psbl_tests+1}='None/Manual';
    
    % Text
    uicontrol(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0.275 0.115 0.2 0.05 ],...
        'String','Current Test Result',...
        'fontsize',12, ...
        'ToolTipString','Set of t-test results to visualize.', ...
        'Style','text');
    % Create sortvar browse button
    dat.h_ptest=uicontrol(dat.fig_id,...
        'CallBack','gui_erp(''change test'');',...
        'Units','normalized', ...
        'Position',[ 0.26 0.075 0.24 0.055 ],...
        'fontsize',12, ...
        'String',ttest_str, ...
        'Value',crnt_ttest, ...
        'Style','popup', ...
        'ToolTipString','Set of t-test results to visualize.', ...
        'tag','bin');

    %%% Time window in which hypotheses tests were done
    %Text
    uicontrol(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0.505 0.115 0.26 0.05 ],...
        'String','Test Window(s) [Min Max]:',...
        'fontsize',12, ...
        'Style','text');
    %Edit Box
    wind_edges=[];
    for nw=1:dat.n_wind,
        if nw==dat.n_wind,
            wind_edges=[wind_edges num2str(test_wind(nw,:))];
        else
            wind_edges=[wind_edges num2str(test_wind(nw,:)) '; ']; 
        end
    end
    dat.h_testwind=uicontrol(dat.fig_id,...
        'CallBack','gui_erp(''update dashed lines'');',...
        'Units','normalized', ...
        'fontsize',12, ...
        'BackGroundColor','w', ...
        'Position',[ 0.525 0.095 0.22 0.04 ],...
        'String',wind_edges, ...
        'Style','edit', ...
        'ToolTipString','Miniumum and maximum time of time window in which hypothesis tests were done.', ...
        'tag','timerange'); 
    
    
    %%% Critical t-score here
    %Text
    uicontrol(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0.76 0.115 0.23 0.05 ],...
        'String','Critical t-score(s)',...
        'fontsize',12, ...
        'Style','text');
    %Edit Box
    dat.h_critval=uicontrol(dat.fig_id,...
        'CallBack','gui_erp(''update critical t'');',...
        'Units','normalized', ...
        'fontsize',12, ...
        'BackGroundColor','w', ...
        'Position',[ 0.79 0.09 0.16 0.045 ],...
        'String',num2str(critical_t), ...
        'Style','edit', ...
        'ToolTipString','Critical t-score(s) for reliable deviation from 0 voltage.', ...
        'tag','timerange');
    
    %%%% WHICH STAT TO PLOT
    stats={'t-score of ERPs','standard error of ERPs','global field power of ERPs'};
    % Text
    uicontrol(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0.01 0.183 0.09 0.03 ],...
        'String','Statistic:',...
        'fontsize',12, ...
        'ToolTipString','Select what to plot in secondary waveform axis: ERP t-scores, standard error of the mean, or global field power.', ...
        'Style','text');
    
    if strcmpi(p.Results.stat,'gfp') || isinf(dat.t_scores(1,1,bin)) || isnan(dat.t_scores(1,1,bin))
        %isinf(dat.t_scores(1,1,bin)) will
        %equal 1 if there are insufficient subjects to calculate t-scores
        ini_val=3;
    elseif strcmpi(p.Results.stat,'t'),
        ini_val=1;
    elseif strcmpi(p.Results.stat,'stder'),
        ini_val=2;
    else
        error('Unrecognized value of input argument ''stat''.');
    end
    dat.h_stat=uicontrol(dat.fig_id,...
        'CallBack','gui_erp(''change stat'');',...
        'Units','normalized', ...
        'Position',[ 0.095 0.165 0.16 0.05 ],...
        'fontsize',12, ...
        'value',ini_val, ...
        'String',stats, ...
        'Style','popup', ...
        'ToolTipString','Plot ERP t-scores or standard error of the mean in secondary axis?', ...
        'tag','stat');
    
    %%% TIME RANGE
    %Text
    uicontrol(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0.01 0.035 0.24 0.05 ],...
        'String','Time Range [Min Max]',...
        'fontsize',12, ...
        'Style','text');
    %Edit Box
    dat.h_timerange=uicontrol(dat.fig_id,...
        'CallBack','gui_erp(''new time limits'');',...
        'Units','normalized', ...
        'fontsize',12, ...
        'BackGroundColor','w', ...
        'Position',[ 0.06 0.013 0.11 0.045 ],...
        'String',num2str(dat.plt_times), ...
        'Style','edit', ...
        'ToolTipString','Miniumum and maximum time on Time x ERP axis.', ...
        'tag','timerange');   
    
    %%% ERP uV RANGE BOX (AXIS A)
    %Text
    uicontrol(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0.26 0.035 0.26 0.05 ],...
        'String','ERP Range [Min Max]',...
        'HorizontalAlignment','left', ...
        'fontsize',12, ...
        'Style','text');
    %Edit Box
    dat.h_statrangeA=uicontrol(dat.fig_id,...
        'CallBack','gui_erp(''new statistic limitsA'');',...
        'Units','normalized', ...
        'fontsize',12, ...
        'BackGroundColor','w', ...
        'Position',[ 0.28 0.013 0.16 0.045 ],...
        'String',num2str(round(vA(3:4)*100)/100), ...
        'ToolTipString','Miniumum and maximum voltage on Time x ERP axis.', ...
        'Style','edit', ...
        'tag','statrange'); 
    
    %%% t-score/StdErr RANGE BOX (AXIS B)
    %Text
    uicontrol(dat.fig_id,...
        'Units','normalized', ...
        'Position',[ 0.50 0.035 0.26 0.05 ],...
        'String','t/StdErr Range [Min Max]',...
        'HorizontalAlignment','left', ...
        'fontsize',12, ...
        'Style','text');
    %Edit Box
    dat.h_statrangeB=uicontrol(dat.fig_id,...
        'CallBack','gui_erp(''new statistic limitsB'');',...
        'Units','normalized', ...
        'fontsize',12, ...
        'BackGroundColor','w', ...
        'Position',[ 0.54 0.013 0.16 0.045 ],...
        'String',num2str(round(vB(3:4)*100)/100), ...
        'ToolTipString','Miniumum and maximum t-scores/voltage on Time x t-Scores/Standard Error/GFP axis.', ...
        'Style','edit', ...
        'tag','statrange');   
    
    % Button to close GUI
    uicontrol(dat.fig_id,...
        'CallBack','close(gcbf);',...
        'Units','normalized', ...
        'Position',[ 0.75 0.012 0.11 0.05 ],...
        'String','Close',...
        'fontsize',12, ...
        'Tag','close', ...
        'backgroundcolor','m', ...
        'ToolTipString','Close GUI', ...
        'Style','pushbutton');
    
    dat.help_msg=sprintf(['Hold mouse cursor over a GUI control for an explanation of what it does.\n\n', ...
        'Time x waveform axis visualizes ERPs, t-score of ERPs, standard error of ERPs, or ERP global field power simultaneously at multiple electrodes.  ' ... 
        'Each colored line corresponds to a different electrode.\n\nClick on Time x uV/t-score/stderr/GFP axis to visualize the scalp topography\n', ...
        'at that point in time.\n\nClick on electrodes in scalp topography to see electrode name.\n\n', ...
        'This GUI was produced by gui_erp.m']);
    
    % Button for help
    uicontrol(dat.fig_id,...
        'CallBack','dat=get(gcbf,''userdata''); helpdlg(dat.help_msg,''ERP GUI Help'');', ...
        'Units','normalized', ...
        'Position',[ 0.87 0.01 0.11 0.05 ],...
        'String','Help',...
        'fontsize',12, ...
        'Tag','help', ...
        'ToolTipString','Click for help', ...
        'Style','pushbutton');
    dat.interrupt=1;
    set(dat.fig_id,'userdata',dat);
elseif strcmpi(cmnd_str,'time jump'),
    new_tpt=find_tpt(str2num(get(dat.h_topo_time,'string')), ...
        dat.times);
    if new_tpt<dat.start_pt,
        errordlg('That time is too early.');
    elseif new_tpt>dat.end_pt,
        errordlg('That time is too late.');
    else
        set(dat.h_topo_time,'string',num2str(dat.times(new_tpt)));
        set(dat.h_lineA,'XData',[1 1]*dat.times(new_tpt));
        set(dat.h_lineB,'XData',[1 1]*dat.times(new_tpt));
        set(dat.fig_id,'userdata',dat);
        redraw_topos(dat,0);
        drawnow;
    end    
elseif strcmpi(cmnd_str,'back one'),
    current_tpt=find_tpt(str2double(get(dat.h_topo_time,'string')), ...
        dat.times);
    new_tpt=current_tpt-1;
    if new_tpt<dat.start_pt,
        errordlg('Can''t go back any further.');
    else
        set(dat.h_topo_time,'string',num2str(dat.times(new_tpt)));
        set(dat.h_lineA,'XData',[1 1]*dat.times(new_tpt));
        set(dat.h_lineB,'XData',[1 1]*dat.times(new_tpt));
        set(dat.fig_id,'userdata',dat);
        redraw_topos(dat,0);
        drawnow;
    end
elseif strcmpi(cmnd_str,'forward one'),
    current_tpt=find_tpt(str2double(get(dat.h_topo_time,'string')), ...
        dat.times);
    new_tpt=current_tpt+1;
    if new_tpt>dat.end_pt,
        errordlg('Can''t go forward any further.');
    else
        set(dat.h_topo_time,'string',num2str(dat.times(new_tpt)));
        set(dat.h_lineA,'XData',[1 1]*dat.times(new_tpt));
        set(dat.h_lineB,'XData',[1 1]*dat.times(new_tpt));
        set(dat.fig_id,'userdata',dat);
        redraw_topos(dat,0);
        drawnow;
    end
elseif strcmpi(cmnd_str,'back'),
    interrupt=0;
    dat.interrupt=0;
    set(dat.h_stop,'enable','on');
    set(dat.h_forward,'enable','off');
    set(dat.h_back,'enable','off');
    set(dat.h_forward1,'enable','off');
    set(dat.h_back1,'enable','off');
    set(dat.fig_id,'userdata',dat);
    drawnow;
    current_tpt=find_tpt(str2double(get(dat.h_topo_time,'string')), ...
        dat.times);
    new_tpt=current_tpt-1;
    loop_ct=0;
    while (new_tpt>=dat.start_pt) && (interrupt~=1),
        set(dat.h_topo_time,'string',num2str(dat.times(new_tpt)));
        set(dat.h_lineA,'XData',[1 1]*dat.times(new_tpt));
        set(dat.h_lineB,'XData',[1 1]*dat.times(new_tpt));
        drawnow('update');
        redraw_topos(dat,0);
        drawnow;
        new_tpt=new_tpt-1;
        loop_ct=loop_ct+1;
        if loop_ct==5,
            loop_ct=0;
            new_tmp=get(dat.fig_id,'userdata');
            interrupt=new_tmp.interrupt;
        end
    end
    drawnow;
    dat.interrupt=0;
    new_tpt=new_tpt+1; %undo step that broke loop
    set(dat.h_topo_time,'string',num2str(dat.times(new_tpt)));
    set(dat.h_stop,'enable','off');
    set(dat.h_forward,'enable','on');
    set(dat.h_back,'enable','on');
    set(dat.h_forward1,'enable','on');
    set(dat.h_back1,'enable','on');
    set(dat.fig_id,'userdata',dat);
elseif strcmpi(cmnd_str,'forward'),
    interrupt=0;
    dat.interrupt=0;
    set(dat.h_stop,'enable','on');
    set(dat.h_forward,'enable','off');
    set(dat.h_back,'enable','off');
    set(dat.h_forward1,'enable','off');
    set(dat.h_back1,'enable','off');
    set(dat.fig_id,'userdata',dat);
    drawnow;
    current_tpt=find_tpt(str2double(get(dat.h_topo_time,'string')), ...
        dat.times);
    new_tpt=current_tpt+1;
    loop_ct=0;
    while (new_tpt<=dat.end_pt) && (interrupt~=1),
        set(dat.h_topo_time,'string',num2str(dat.times(new_tpt)));
        set(dat.h_lineA,'XData',[1 1]*dat.times(new_tpt));
        set(dat.h_lineB,'XData',[1 1]*dat.times(new_tpt));
        drawnow('update');
        redraw_topos(dat,0);
        drawnow;
        new_tpt=new_tpt+1;
        loop_ct=loop_ct+1;
        if loop_ct==5,
            loop_ct=0;
            new_tmp=get(dat.fig_id,'userdata');
            interrupt=new_tmp.interrupt;
        end
    end
    drawnow;
    dat.interrupt=0;
    new_tpt=new_tpt-1; %undo step that broke loop
    set(dat.h_topo_time,'string',num2str(dat.times(new_tpt)));
    set(dat.h_stop,'enable','off');
    set(dat.h_forward,'enable','on');
    set(dat.h_back,'enable','on');
    set(dat.h_forward1,'enable','on');
    set(dat.h_back1,'enable','on');
    set(dat.fig_id,'userdata',dat);
elseif strcmpi(cmnd_str,'change stat'),
    stat=get(dat.h_stat,'value');
    current_bin=get(dat.h_bin,'value');
    if (stat~=3),
        %make sure there are sufficients subjects for that bin
        if isinf(dat.t_scores(1,1,current_bin)) || isnan(dat.t_scores(1,1,current_bin)) 
            stat=3;
            set(dat.h_stat,'value',stat);
            if (stat==2)
                warndlg('Only one participant contributed to this bin. Thus the grand average standard error cannot be esimated.','gui_erp Warning');
            else
                warndlg('Only one participant contributed to this bin. Thus the grand average t-scores cannot be derived.','gui_erp Warning');
            end
        end
    end
        
    if (stat==1),
        % t-scores
        if ~isempty(dat.h_crit1),
            delete(dat.h_crit1);
        end
        if ~isempty(dat.h_crit2),
            delete(dat.h_crit2);
        end
        dat.h_crit1=[];
        dat.h_crit2=[];
        
        if ~isempty(dat.h_alph1),
            delete(dat.h_alph1);
        end
        if ~isempty(dat.h_alph2),
            delete(dat.h_alph2);
        end
        dat.h_alph1=[];
        dat.h_alph2=[];
    end
    
    %redraw time plot
    draw_waveforms(dat,0);
    gui_erp('update dashed lines'); %also redraws topo  
elseif strcmpi(cmnd_str,'change bin'),
    new_bin=get(dat.h_bin,'value');
    crnt_ttest_id=get(dat.h_ptest,'value');
    n_t_tests=length(dat.t_tests);
    if (crnt_ttest_id==(n_t_tests+1)) || (dat.t_tests(crnt_ttest_id).bin~=new_bin),
        % If currently selected set of t-tests is for a different bin
        % than is currently visualized search for a t-test set on
        % this bin. Use the first test if more than one test has
        % been performed on the bin.  If no test is
        % found, set to None/Manual.
        crnt_ttest_id=n_t_tests+1;
        for a=1:n_t_tests,
            if dat.t_tests(a).bin==new_bin;
                crnt_ttest_id=a;
                break;
            end
        end
        set(dat.h_ptest,'value',crnt_ttest_id);
        dat=update_test(dat);
    end
    draw_waveforms(dat,1);
    gui_erp('update dashed lines'); %also redraws topo plot
elseif strcmpi(cmnd_str,'new time limits'),
    %time limits
    plotted_times=dat.times(dat.start_pt:dat.end_pt);
    tme_lim=str2num(get(dat.h_timerange,'string'));
    if length(tme_lim)~=2,
        errordlg('Enter exactly two values for time range (desired Min and Max values).','gui_erp Error');
        set(dat.h_timerange,'string',num2str([plotted_times(1) plotted_times(end)])); %reset time limits
    elseif tme_lim(2)<tme_lim(1),
        errordlg('The second time range value must be greater than the first (enter desired Min then Max).','gui_erp Error');
        set(dat.h_timerange,'string',num2str([plotted_times(1) plotted_times(end)])); %reset time limits
    elseif (tme_lim(1)<dat.times(1)) || (tme_lim(2)>dat.times(end))
        errordlg(sprintf('Time range values must be between %d and %d msec.',dat.times(1),dat.times(end)),'gui_erp Error');
        set(dat.h_timerange,'string',num2str([plotted_times(1) plotted_times(end)])); %reset time limits
    else
        %change time limits
        set(dat.h_timeA,'xlim',tme_lim);
        set(dat.h_timeB,'xlim',tme_lim);
        dat.start_pt=find_tpt(tme_lim(1),dat.times);
        dat.end_pt=find_tpt(tme_lim(2),dat.times);
        
        %update max abs val
        dat.absmxA=max(max(abs(dat.showingA(:,dat.start_pt:dat.end_pt))));
        dat.absmxB=max(max(abs(dat.showingB(:,dat.start_pt:dat.end_pt))));
        current_tpt=find_tpt(str2double(get(dat.h_topo_time,'string')), ...
            dat.times);
        if (current_tpt>dat.end_pt) || (current_tpt<dat.start_pt)
           % Time point whose topography is currently being visualized
           % falls outside of time range.  Reset current time point to peak
           % t-score latency within visualized time window.
           [dummy, mx_tpt]=max(max(abs(dat.showing_t(:,dat.start_pt:dat.end_pt))));
           showing_pts=dat.start_pt:dat.end_pt;
           set(dat.h_lineA,'xdata',[1 1]*dat.times(showing_pts(mx_tpt)));
           set(dat.h_lineB,'xdata',[1 1]*dat.times(showing_pts(mx_tpt)));
           
           set(dat.h_topo_time,'string',num2str(dat.times(showing_pts(mx_tpt))));
        end
        
        redraw_topos(dat,1);
        set(dat.fig_id,'userdata',dat);
    end
elseif strcmpi(cmnd_str,'new statistic limitsA'),
    %statistic limits
    plotted_limits=get(dat.h_timeA,'ylim');
    new_lim=str2num(get(dat.h_statrangeA,'string'));
    if length(new_lim)~=2,
        errordlg('Enter exactly two values for statistic range (desired Min and Max values).','gui_erp Error');
        set(dat.h_statrangeA,'string',num2str([plotted_limits(1) plotted_limits(end)])); %reset statistic limits
    elseif new_lim(2)<new_lim(1),
        errordlg('The second statistic range value must be greater than the first (enter desired Min then Max).','gui_erp Error');
        set(dat.h_statrangeA,'string',num2str([plotted_limits(1) plotted_limits(end)])); %reset statistic limits
    else
        %change time line limits
        set(dat.h_timeA,'ylim',new_lim);
        set(dat.h_t0lineA,'YData',new_lim); %adjust time=0 line
        set(dat.h_lineA,'YData',new_lim); %adjust current time pt line
        
        %change dashed time line limits (hypothesis test window)
        if ~isempty(dat.h_wind1A)
           set(dat.h_wind1A,'YData',new_lim); 
           set(dat.h_wind2A,'YData',new_lim);
        end
    end
elseif strcmpi(cmnd_str,'new statistic limitsB'),
    %statistic limits
    plotted_limits=get(dat.h_timeB,'ylim');
    new_lim=str2num(get(dat.h_statrangeB,'string'));
    if length(new_lim)~=2,
        errordlg('Enter exactly two values for statistic range (desired Min and Max values).','gui_erp Error');
        set(dat.h_statrangeB,'string',num2str([plotted_limits(1) plotted_limits(end)])); %reset statistic limits
    elseif new_lim(2)<new_lim(1),
        errordlg('The second statistic range value must be greater than the first (enter desired Min then Max).','gui_erp Error');
        set(dat.h_statrangeB,'string',num2str([plotted_limits(1) plotted_limits(end)])); %reset statistic limits
    else
        %change time line limits
        set(dat.h_timeB,'ylim',new_lim);
        set(dat.h_t0lineB,'YData',new_lim); %adjust time=0 line
        set(dat.h_lineB,'YData',new_lim); %adjust current time pt line
        
        %change dashed time line limits (hypothesis test window)
        if ~isempty(dat.h_wind1B)
            set(dat.h_wind1B,'YData',new_lim);
            set(dat.h_wind2B,'YData',new_lim);
        end
    end
elseif strcmpi(cmnd_str,'redraw topo'),
    redraw_topos(dat,1); %we need this command string for buttdownfcn calls
elseif strcmpi(cmnd_str,'update dashed lines'),
    test_wind=str2num(get(dat.h_testwind,'string'));
    
    if (size(test_wind,2)>2) || (isempty(test_wind) && ~isempty(get(dat.h_testwind,'string')))
        errordlg('Enter pairs of values for time window range (desired Min and Max values) or leave box blank. Separate multiple pairs with semicolons.', ...
            'gui_erp Error');
    elseif length(dat.critical_t)>2,
        errordlg('Enter one (for one-tailed test) or two (for two-tailed test) value(s) for critical t-scores or leave box blank.', ...
            'gui_erp Error');
    else
        if isempty(test_wind),
            test_wind=[NaN NaN];
        end
        
        % remove existing critical value markers & test time window markers
        if ~isempty(dat.h_crit1),
            delete(dat.h_crit1);
        end
        if ~isempty(dat.h_crit2),
            delete(dat.h_crit2);
        end
        if ~isempty(dat.h_alph1),
            delete(dat.h_alph1);
        end
        if ~isempty(dat.h_alph2),
            delete(dat.h_alph2);
        end
        if ~isempty(dat.h_wind1A),
            delete(dat.h_wind1A);
        end
        if ~isempty(dat.h_wind2A),
            delete(dat.h_wind2A);
        end
        if ~isempty(dat.h_wind1B),
            delete(dat.h_wind1B);
        end
        if ~isempty(dat.h_wind2B),
            delete(dat.h_wind2B);
        end
        dat.h_wind1A=[];
        dat.h_wind2A=[];
        dat.h_wind1B=[];
        dat.h_wind2B=[];
        dat.h_crit1=[];
        dat.h_crit2=[];
        dat.h_alph1=[];
        dat.h_alph2=[];
        dat.n_wind=size(test_wind,1);
        
        if ~sum(isnan(test_wind))
            y_rngA=get(dat.h_timeA,'YLim');
            y_rngB=get(dat.h_timeB,'YLim');
            
            axes(dat.h_timeA); %warning message says that I should give axes an output to speed it up, but axes doesn't produce an output
            for nw=1:dat.n_wind,
                dat.h_wind1A(nw)=plot([1 1]*test_wind(nw,1),y_rngA,'k--');
                set(dat.h_wind1A(nw),'linewidth',2);
                dat.h_wind2A(nw)=plot([1 1]*test_wind(nw,2),y_rngA,'k--');
                set(dat.h_wind2A(nw),'linewidth',2);
            end
            
            axes(dat.h_timeB);
            for nw=1:dat.n_wind,
                dat.h_wind1B(nw)=plot([1 1]*test_wind(nw,1),y_rngB,'k--');
                set(dat.h_wind1B(nw),'linewidth',2);
                dat.h_wind2B(nw)=plot([1 1]*test_wind(nw,2),y_rngB,'k--');
                set(dat.h_wind2B(nw),'linewidth',2);
            end
            
            %Find out if visualized units are microvolts or t-scores
            stat=get(dat.h_stat,'value');
            
            % if stat=1, t-scores, else data are in microvolts (for stderr)
            if ~sum(isnan(dat.critical_t)) && (stat==1) && ~isempty(dat.critical_t)
                for nw=1:dat.n_wind,
                    dat.h_crit1(nw)=plot(test_wind(nw,:),[1 1]*dat.critical_t(1),'r--');
                    set(dat.h_crit1(nw),'linewidth',3);
                    if nw==1,
                        %write alpha value on figure
                        tm_rng=dat.times(end)-dat.times(1);
                        if strcmpi(dat.mltplcmp_crct,'fdr')
                            dat.h_alph1=text(test_wind(nw,1)-tm_rng*.02,dat.critical_t(1), ...
                                ['q=' num2str(rnd_orderofmag(dat.alpha))]);
                        else
                            dat.h_alph1=text(test_wind(nw,1)-tm_rng*.02,dat.critical_t(1), ...
                                ['\alpha=' num2str(rnd_orderofmag(dat.alpha))]);
                        end
                        set(dat.h_alph1,'color','r','fontweight','normal','fontsize',12, ...
                            'horizontalalignment','right','backgroundcolor',[1 1 1], ...
                            'edgecolor',[1 1 1]*.3,'clipping','on','fontname','fixedwidth');
                    end
                    if length(dat.critical_t)>1,
                        dat.h_crit2(nw)=plot(test_wind(nw,:),[1 1]*dat.critical_t(2),'r--');
                        set(dat.h_crit2(nw),'linewidth',3);
                        if nw==1,
                            %write alpha value on figure
                            if strcmpi(dat.mltplcmp_crct,'fdr')
                                dat.h_alph2=text(test_wind(1)-tm_rng*.02,dat.critical_t(2), ...
                                    ['q=' num2str(rnd_orderofmag(dat.alpha))]);
                            else
                                dat.h_alph2=text(test_wind(1)-tm_rng*.02,dat.critical_t(2), ...
                                    ['\alpha=' num2str(rnd_orderofmag(dat.alpha))]);
                            end
                            set(dat.h_alph2,'color','r','fontweight','normal','fontsize',12, ...
                                'horizontalalignment','right','backgroundcolor',[1 1 1], ...
                                'edgecolor',[1 1 1]*.3,'clipping','on','fontname','fixedwidth');
                        end
                    end
                end
            end
        end
        set(dat.fig_id,'userdata',dat);
        redraw_topos(dat,1);
    end
elseif strcmpi(cmnd_str,'update critical t'),
    critical_t=str2num(get(dat.h_critval,'string'));
    if ~isequal(sort(critical_t),sort(dat.critical_t)),
        new_alpha=[];
        while isempty(new_alpha),
            new_alpha=inputdlg({'Enter alpha level for new critical t-score(s)'}, ...
                'Set New Alpha Level',1,{num2str(dat.alpha)});
            new_alpha=str2num(new_alpha{1});
            if isempty(new_alpha),
                errrodlg('You must enter a numeric value for the new alpha level','gui_erp Error');
            end
        end
        dat.alpha=new_alpha;
        dat.critical_t=critical_t;
    end
    set(dat.fig_id,'userdata',dat);
    gui_erp('update dashed lines');
elseif strcmpi(cmnd_str,'change test'),
    [dat, redraw_wform]=update_test(dat);    
    if redraw_wform,
        draw_waveforms(dat,1);
    end
    gui_erp('update dashed lines'); %also redraws topo
else
    errordlg(sprintf('Command "%s" not recognized by gui_erp.m',cmnd_str), ...
        'gui_erp Error');
end


function draw_waveforms(dat,redrawA)
% Redraws the channel waveforms (i.e., ERPs, t-scores of ERPs, standard error of ERPs, or GFP of ERPs)
% and the color bar the topography.
%
% Used by: 'change bin', 'change stat', and 'change test'
%
% Input:
%   dat     - The userdata of the GUI
%   redrawA - If non-zero, axis A (the ERP waveformss) will also be redrawn

new_bin=get(dat.h_bin,'value');
stat=get(dat.h_stat,'value');
crnt_ttest_id=get(dat.h_ptest,'value');
n_t_tests=length(dat.t_tests);
if crnt_ttest_id==n_t_tests+1,
    %no t-test results visualized
    dat.showing_chans=dat.loaded_chans;
else
    dat.showing_chans=intersect(dat.t_tests(crnt_ttest_id).used_chan_ids, ...
        dat.loaded_chans);
end

%redraw Axis B time plot
axes(dat.h_timeB);
delete(dat.h_showingB);
delete(dat.h_lineB);
plotted_times=dat.times(dat.start_pt:dat.end_pt);
%store t-scores for waveform that will be shown
if (crnt_ttest_id<=n_t_tests) && dat.t_tests(crnt_ttest_id).null_mean
    %The mean of the null hypothesis is non-zero
    dat.showing_t=squeeze( (dat.erp(dat.showing_chans,:,new_bin)- ...
        dat.t_tests(crnt_ttest_id).null_mean)./dat.stder(dat.showing_chans,:,new_bin) );
else
    dat.showing_t=squeeze(dat.t_scores(dat.showing_chans,:,new_bin));
end

if (stat==1),
    dat.showingB=dat.showing_t;
    cbar_title='t-score';
    ylab='t-score';
elseif (stat==2),
    dat.showingB=squeeze(dat.stder(dat.showing_chans,:,new_bin));
    cbar_title='\muV';
    ylab='\muV (StdEr)';
elseif (stat==3),
    dat.showingB=dat.gfp(:,new_bin)';
    cbar_title='\muV';
    ylab='\muV (GFP)';
end
if size(dat.showingB,1)==1,
    %only one channel being plot
    dat.h_showingB=plot(dat.times,dat.showingB);
    cbar_title='Not Applicable';
    cbar_title_fontsize=10;
else
    dat.h_showingB=plot(dat.times,dat.showingB');
    cbar_title_fontsize=14;
end
set(dat.h_time_ylabB,'string',ylab);

stat_mx=max(max(dat.showingB));
stat_mn=min(min(dat.showingB));
%Make sure new bin has data in it
if isnan(stat_mx) || isnan(stat_mn)
   error('Bin %d does not appear to have any data in it.',new_bin); 
end
stat_rng=stat_mx-stat_mn;
stat_plt_rng=[stat_mn-stat_rng*.02 stat_mx+stat_rng*.02];
stat_plt_rng=round(stat_plt_rng*100)/100;
axis([plotted_times(1) plotted_times(end) stat_plt_rng]);
set(dat.h_statrangeB,'string',num2str(stat_plt_rng));
set(dat.h_t0lineB,'YData',stat_plt_rng); %adjust time=0 line
current_tpt=find_tpt(str2double(get(dat.h_topo_time,'string')), ...
    dat.times);
h_lineB=plot([1 1]*dat.times(current_tpt),stat_plt_rng,'k');
set(h_lineB,'linewidth',2);
dat.h_lineB=h_lineB;

bdf_code = [ 'tmppos = get(gca, ''currentpoint'');' ...
    'dat=get(gcbf, ''userdata'');' ...
    'new_tpt=find_tpt(tmppos(1,1),dat.times);' ...
    'set(dat.h_lineA,''XData'',[1 1]*dat.times(new_tpt));' ...
    'set(dat.h_lineB,''XData'',[1 1]*dat.times(new_tpt));' ...
    'set(dat.h_topo_time,''string'',num2str(dat.times(new_tpt)));' ...
    'set(dat.fig_id,''userdata'',dat);' ...
    'gui_erp(''redraw topo'');' ...
    'drawnow;' ...
    'clear latpoint dattmp tmppos;' ...
    ];
set(dat.h_showingB,'ButtonDownFcn',bdf_code);
set(dat.h_timeB,'ButtonDownFcn',bdf_code);

%redraw topo color bar for Axis B
dat.absmxB=max(max(abs(dat.showingB(:,dat.start_pt:dat.end_pt))));
axes(dat.h_cbarB);
cla;
cbar(dat.h_cbarB);
absmx=round(dat.absmxB*100)/100;
set(gca,'xticklabel',[-absmx 0 absmx]);
h_cbar_title=title(cbar_title);
set(h_cbar_title,'fontsize',cbar_title_fontsize);

if redrawA,
    %redraw Axis A time plot
    new_title=['Bin ' int2str(new_bin) ': ' dat.bindesc{new_bin}];
    title_max_char=43;
    if length(new_title)>title_max_char,
        new_title=new_title(1:title_max_char);
    end
    set(dat.h_time_title,'string',new_title);
    
    
    %redraw Axis A time plot
    axes(dat.h_timeA);
    delete(dat.h_showingA);
    delete(dat.h_lineA);
    plotted_times=dat.times(dat.start_pt:dat.end_pt);
    %store t-scores for waveform that will be shown
    if (crnt_ttest_id<=n_t_tests) && dat.t_tests(crnt_ttest_id).null_mean
        %The mean of the null hypothesis is non-zero
        dat.showing_t=squeeze( (dat.erp(dat.showing_chans,:,new_bin)- ...
            dat.t_tests(crnt_ttest_id).null_mean)./dat.stder(dat.showing_chans,:,new_bin) );
    else
        dat.showing_t=squeeze(dat.t_scores(dat.showing_chans,:,new_bin));
    end
    
    dat.showingA=squeeze(dat.erp(dat.showing_chans,:,new_bin));
    cbar_title='\muV';
    ylab='\muV (ERP)';
    if size(dat.showingA,1)==1,
        %only one channel being plot
        dat.h_showingA=plot(dat.times,dat.showingA);
        cbar_title='Not Applicable';
        cbar_title_fontsize=10;
    else
        dat.h_showingA=plot(dat.times,dat.showingA');
    end
    set(dat.h_time_ylabA,'string',ylab);
    
    stat_mx=max(max(dat.showingA));
    stat_mn=min(min(dat.showingA));
    stat_rng=stat_mx-stat_mn;
    stat_plt_rng=[stat_mn-stat_rng*.02 stat_mx+stat_rng*.02];
    stat_plt_rng=round(stat_plt_rng*100)/100;
    axis([plotted_times(1) plotted_times(end) stat_plt_rng]);
    set(dat.h_statrangeA,'string',num2str(stat_plt_rng));
    set(dat.h_t0lineA,'YData',stat_plt_rng); %adjust time=0 line
    current_tpt=find_tpt(str2double(get(dat.h_topo_time,'string')), ...
        dat.times);
    h_lineA=plot([1 1]*dat.times(current_tpt),stat_plt_rng,'k');
    set(h_lineA,'linewidth',2);
    dat.h_lineA=h_lineA;
    
    bdf_code = [ 'tmppos = get(gca, ''currentpoint'');' ...
        'dat=get(gcbf, ''userdata'');' ...
        'new_tpt=find_tpt(tmppos(1,1),dat.times);' ...
        'set(dat.h_lineA,''XData'',[1 1]*dat.times(new_tpt));' ...
        'set(dat.h_lineB,''XData'',[1 1]*dat.times(new_tpt));' ...
        'set(dat.h_topo_time,''string'',num2str(dat.times(new_tpt)));' ...
        'set(dat.fig_id,''userdata'',dat);' ...
        'gui_erp(''redraw topo'');' ...
        'drawnow;' ...
        'clear latpoint dattmp tmppos;' ...
        ];
    set(dat.h_timeA,'ButtonDownFcn',bdf_code);
    set(dat.h_showingA,'ButtonDownFcn',bdf_code);
    
    %redraw topo color bar for Axis A
    dat.absmxA=max(max(abs(dat.showingA(:,dat.start_pt:dat.end_pt))));
    axes(dat.h_cbarA);
    cla;
    cbar(dat.h_cbarA);
    absmx=round(dat.absmxA*100)/100;
    set(gca,'xticklabel',[-absmx 0 absmx]);
    h_cbar_title=title(cbar_title);
    set(h_cbar_title,'fontsize',cbar_title_fontsize);
end

set(dat.fig_id,'userdata',dat);

%note, redraw_topos is called by update dashed lines



function [dat, redraw_wform]=update_test(dat)
% Updates the t-test results paramters (e.g., critical t-score)
% Used by: 'change test' & 'change bin'
%
% Inputs:
%   dat - The userdata of the GUI
%
% Outputs:
%   dat            - The revised of the GUI
%   redraw_wform   - 1 if the waveforms in the GUI need to be redrawn
%                    (e.g., if the bin being visualized has changed),
%                    otherwise it's 0
%

crnt_ttest=get(dat.h_ptest,'value');

n_psbl_tests=length(dat.psbl_tests);
crnt_bin=get(dat.h_bin,'value');
redraw_wform=0;
if n_psbl_tests<crnt_ttest,
    %manual/no test, remove all test info/visualizations
    set(dat.h_testwind,'string',[]); %get rid of test wind
    dat.alpha=[];
    dat.critical_t=[];
    dat.mltplcmp_crct=[];
    set(dat.h_critval,'string',[]);
    
    if ~isequal(dat.loaded_chans,dat.showing_chans)
        redraw_wform=1;
    end    
else
    ptest_bin=dat.t_tests(crnt_ttest).bin;
    if ptest_bin~=crnt_bin,
        redraw_wform=1;
        set(dat.h_bin,'value',ptest_bin);
    end
    
    %If the mean of the null hypothesis of the previously shown test
    %results differs from that of the new test results, redraw the waveform
    crnt_null_mean=dat.null_mean;
    if crnt_null_mean~=dat.t_tests(crnt_ttest).null_mean,
        redraw_wform=1;
    end
    dat.null_mean=dat.t_tests(crnt_ttest).null_mean;
    
    %test time window(s)
    test_wind=dat.t_tests(crnt_ttest).time_wind;
    dat.n_wind=size(test_wind,1);
    wind_edges=[];
    for nw=1:dat.n_wind,
        if nw==dat.n_wind,
            wind_edges=[wind_edges num2str(test_wind(nw,:))];
        else
            wind_edges=[wind_edges num2str(test_wind(nw,:)) '; '];
        end
    end
    set(dat.h_testwind,'string',wind_edges);
    
    % if test was performed including channels that are not currently
    % loaded, throw a warning
    if ~isempty(setdiff(dat.t_tests(crnt_ttest).used_chan_ids, ...
            dat.loaded_chans)),
        warndlg('The set of t-tests you selected includes some channels that were excluded when you created this GUI.  Keep this in mind when interpreting the visualization.','gui_erp Warning');
    end
    
    % if test was performed on a set of channels that differs from what
    % is currently shown, replot the data with the 'change bin' command
    if ~isequal(dat.t_tests(crnt_ttest).used_chan_ids,dat.showing_chans)
        redraw_wform=1;
    end
    
    % set critical values and alpha level
    if isnan(dat.t_tests(crnt_ttest).estimated_alpha)
        dat.alpha=dat.t_tests(crnt_ttest).desired_alphaORq;
        dat.mltplcmp_crct='fdr';
    else
        dat.alpha=dat.t_tests(crnt_ttest).estimated_alpha;
        dat.mltplcmp_crct='perm';
    end
    dat.critical_t=dat.t_tests(crnt_ttest).crit_t;
    set(dat.h_critval,'string',num2str(dat.critical_t));
end
set(dat.fig_id,'userdata',dat);


function redraw_topos(dat,redraw_cbars)

current_tpt=find_tpt(str2double(get(dat.h_topo_time,'string')), ...
    dat.times);

test_wind=str2num(get(dat.h_testwind,'string'));
critical_t=str2num(get(dat.h_critval,'string'));

if ~isempty(test_wind) && ~isempty(critical_t),
    %if current time point is in a test window look for channels with
    %sig effects
    crnt_ttest=get(dat.h_ptest,'value');
    current_tpt_in_ms=dat.times(current_tpt);
    n_showing=length(dat.showing_chans);
    sig_chans=zeros(1,n_showing);    
    for nw=1:dat.n_wind,
        if (current_tpt_in_ms>=test_wind(nw,1)) && (current_tpt_in_ms<=test_wind(nw,2)),
            %if isnan(critical_t) YODA
            if ~strcmpi(dat.mltplcmp_crct,'fdr') && (length(critical_t)==1) && isnan(critical_t) % critical_t is always the scalar NaN for cluster based permutation tests
                %cluster based test
                pval_tpt_id=find(dat.t_tests(crnt_ttest).used_tpt_ids==current_tpt);
                if ~isempty(pval_tpt_id)
                    %if there's a significant electrode
                    sig_chans_temp=find(dat.t_tests(crnt_ttest).adj_pval(:,pval_tpt_id)<dat.t_tests(crnt_ttest).desired_alphaORq);
                    sig_chans_temp=dat.t_tests(crnt_ttest).used_chan_ids(sig_chans_temp); %convert sig channel indices into channel indices in the original GND/GRP variable
                    sig_chans=zeros(1,n_showing);
                    for a=1:n_showing,
                        if ismember(dat.showing_chans(a),sig_chans_temp)
                            sig_chans(a)=1;
                        end
                    end
                end
            else
                if length(critical_t)==2,
                    % I don't think one critical t value can be NaN and not
                    % the other, but just in case
                    if isnan(critical_t(2))
                        sig_chans=[];
                    else
                        sig_chans(dat.showing_t(:,current_tpt)>max(critical_t))=1;
                    end
                    if ~isnan(critical_t(1))
                        sig_chans(dat.showing_t(:,current_tpt)<min(critical_t))=1;
                    end
                else
                    if isnan(critical_t) %FDR control con produce NaN critical_t's if not test are significant
                        sig_chans=[];
                    elseif critical_t>0,
                        sig_chans(dat.showing_t(:,current_tpt)>critical_t)=1;
                    else
                        sig_chans(dat.showing_t(:,current_tpt)<critical_t)=1;
                    end
                end
            end
            break; %Visualized time point is in this window. Break out of for loop since there's no need to look at additional time windows
        end
    end
    sig_chans=find(sig_chans>0);
else
    sig_chans=[];
end

% ERP TOPO
set(dat.fig_id,'CurrentAxes',dat.h_topoA); %<-changes current axis without updating GUI graphics
%which is why I use this instead of axes(dat.h_topoA)
if size(dat.showingA,1)<=2,
    %two or fewer channels, we can only plot electrode locations
    topoplotMK(dat.showingA(:,current_tpt),dat.chanlocs(dat.showing_chans), ...
        'style','blank','plain_blank',1,'emarker2',{sig_chans,'o',[1 1 1],4});
    
    if get(dat.h_stat,'value')==3,
        %GFP topo
        if redraw_cbars,
            %statistic has been changed, need to erase old topo
            set(dat.fig_id,'CurrentAxes',dat.h_topoB);
            topoplotMK([],dat.chanlocs(dat.showing_chans), ...
                'style','blank','plain_blank',1);
        end
    else
        % t-score/stderr
        set(dat.fig_id,'CurrentAxes',dat.h_topoB);
        topoplotMK(dat.showingB(:,current_tpt),dat.chanlocs(dat.showing_chans), ...
            'style','blank','plain_blank',1,'emarker2',{sig_chans,'o',[1 1 1],4});
    end
else
    topoplotMK(dat.showingA(:,current_tpt),dat.chanlocs(dat.showing_chans),'maplimits', ...
        [-1 1]*dat.absmxA,'emarker2',{sig_chans,'o',[1 1 1],4});
    set(findobj(gca,'type','patch'),'facecolor',[1 1 1]*.702); %make topoplot background color match that of GUI
    
    if get(dat.h_stat,'value')==3,
        %GFP topo
        if redraw_cbars,
            %statistic has been changed, need to erase old topo
            set(dat.fig_id,'CurrentAxes',dat.h_topoB);
            topoplotMK([],dat.chanlocs(dat.showing_chans), ...
                'style','blank','plain_blank',1);
            set(findobj(gca,'type','patch'),'facecolor',[1 1 1]*.702); %make topoplot background color match that of GUI
        end
    else
        % t-score/stderr
        set(dat.fig_id,'CurrentAxes',dat.h_topoB);
        topoplotMK(dat.showingB(:,current_tpt),dat.chanlocs(dat.showing_chans),'maplimits', ...
            [-1 1]*dat.absmxB,'emarker2',{sig_chans,'o',[1 1 1],4});
        set(findobj(gca,'type','patch'),'facecolor',[1 1 1]*.702); %make topoplot background color match that of GUI
    end
end

if (redraw_cbars)
    % ERP COLOR BAR TICK LABELS
    %axes(dat.h_cbarA);
    set(dat.fig_id,'CurrentAxes',dat.h_cbarA);
    absmx=round(dat.absmxA*100)/100;
    set(gca,'xticklabel',[-absmx 0 absmx]);
    
    % t-Score/StdErr/GFP COLOR BAR TICK LABELS
    %axes(dat.h_cbarB);
    set(dat.fig_id,'CurrentAxes',dat.h_cbarB);
    absmx=round(dat.absmxB*100)/100;
    set(gca,'xticklabel',[-absmx 0 absmx]);
end
drawnow

