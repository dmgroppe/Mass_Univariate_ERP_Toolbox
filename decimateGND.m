% decimateGND() - Downsamples a set of ERPs to a lower sampling rate.
%                 Useful for reducing the number of dependent variables and 
%                 conserving memory.
% Usage:
%  >> GND=decimateGND(GND,decfactor,method,bsln_wind,verblevel);
%
% Required Inputs:
%   GND       - A GND structure variable.  To create a GND variable from
%               Kutaslab ERP files (e.g., *.mas files) use avgs2GND.m.  To 
%               do the same from EEGLAB *.set files use sets2GND.m.  See Mass
%               Univariate ERP Toolbox documentation for detailed information  
%               about the format of a GND variable. 
%   decfactor - [positive integer] The factor by which to reduce the data
%               set.  For example, if decfactor is 2, the data will have 
%               approximately half as many time points.
%
% Optional Inputs:
%   method     - ['boxcar' or 'fir'] If 'boxcar,' data are decimated with a
%                boxcar moving average (this is the algorithm used by 
%                Kutaslab's UNIX program "avg").  If 'fir,' data are decimated
%                using the MATLAB function decimate.m and a 30th order FIR
%                filter. See comments in this file (decimateGND.m) for more 
%                information on how exactly decimation is done.  Note, you
%                MUST HAVE THE MATLAB SIGNAL PROCESSING TOOLBOX to use the 
%                'fir' option.  The 'boxcar' option does NOT require the 
%                signal processing toolbox {default: 'boxcar'}
%   bsln_wind  - [vector] Two element vector specifying the beginning and
%                end (in ms) of the baseline time window (e.g., [-100 -4]).
%                The mean amplitude across all time points within and 
%                including those times will be removed from each ERP.  Data
%                should be re-baselined after decimation, because baseline
%                mean amplitude may no longer be zero after decimation.
%                {default: all time points before 0}
%   save_GND   - ['yes' or 'no'] If 'yes', the GND variable will be
%                saved to disk after the permutation test is completed 
%                and added to it. User will first be prompted to verify 
%                file name and path. {default: 'yes'}
%   verblevel  - An integer specifiying the amount of information you want
%                this function to provide about what it is doing during runtime.
%                  Options are:
%                    0 - quiet, only show errors, warnings, and EEGLAB reports
%                    1 - stuff anyone should probably know
%                    2 - stuff you should know the first time you start working
%                        with a data set {default value}
%                    3 - stuff that might help you debug (show all
%                        reports)
%
% Example:
% To convert a GND variable from 250 Hz to 125 Hz using a boxcar window:
% >>GND=decimateGND(GND,2);
%
% To convert a GND variable from 250 Hz to 125 Hz using a FIR anti-aliasing 
% filter and then re-baselining the data:
% >>GND=decimateGND(GND,2,'fir',[-100 0]);
%
%
% Author:
% David Groppe ('boxcar' algorithm from Paul Krewski)
% Kutaslab, 3/2010

function GND=decimateGND(GND,decfactor,method,bsln_wind,save_GND,verblevel)

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 3/16/2011-Function now works when GND.cals field is empty.  This field is
% unique to Kutaslab data and is empty for data obtained from other labs.
% Thanks to Andrew Hill for finding this bug.
%
% 3/23/2011-save_GND option added. 
%
% 3/12/2013-'fir' option now works when GND.cals field is empty (bug fix on
% 3/16/2011 was accidentally incomplete).  Thanks to Aaron Newman for
% reporting this.

%Notes about decimation:
% -When using the 'boxcar' method, the size of the boxcar moving window is
% decfactor+1.  For example, if decfactor=2, each time point is replaced by
% the mean of it and the immediately preceding and following time point.
% Time points at the edge of the ERPs are discarded since they lack
% sufficient preceding or following time points.  Every other time point is
% discarded.  This method has the advantage that the latency of any effect
% won't be distorted (i.e., data at a pre-decimated time point only effects
% the data at a single post-decimated time point).  The DISADVANTAGE of this
% method is that moving average windows have ripply stop bands.  Thus it
% won't do a great job of suppressing high frequency activity (e.g., 60 Hz
% and will surely result in some aliasing.  Unless you have a lot of
% 60 Hz or EMG in your data this might not be a significant problem since 
% ERPs typically have little power at those high frequencies.
%   I suspect that the 'fir' option could potentially shift effects a bit
% in time, but I haven't verified this yet.  When I tried running the
% filter on an impulse (a waveform with a single non-zero value), it showed
% no evidence of this.  In contrast, the Chebyshev filter (decimate.m's
% default filter) does spread an impulse forward and backward in time and
% could distort the latency of effects.
%

if nargin<3,
    method='boxcar';
else
   if ~(strcmpi(method,'boxcar') || strcmpi(method,'fir')),
       error('Arugment ''method'' needs to be ''boxcar'' or ''fir''');
   end
end

if nargin<4
    bsln_wind=[];
elseif ~isempty(bsln_wind),
   if length(bsln_wind)~=2,
       error('Argument bsln_wind needs to be a two element vector.');
   end
end

if nargin<5
    save_GND='yes';
elseif ~(strcmpi(save_GND,'yes') || strcmpi(save_GND,'no'))
    error('Argument save_GND needs to be ''yes'' or ''no''.');
end

global VERBLEVEL;
if nargin<6
    if isempty(VERBLEVEL),
        VERBLEVEL=2;
    end
else
    VERBLEVEL=verblevel;
end


%Erase any t-test results as they won't be valid anymore
if ~isempty(GND.t_tests),
    %but ask first
    if VERBLEVEL>1,
        resp=[];
        while ~strcmpi(resp,'y') && ~strcmpi(resp,'n') && ~strcmpi(resp,'yes') ...
                && ~strcmpi(resp,'no'),
            fprintf('These data have t-test results stored with them that won''t be accurate after the data have been decimated.\nFor this reason, they will be erased.\n');
            resp=input(sprintf('Continue with decimation (t-tests will be erased)? [y or n] '),'s');
        end
        if strcmpi('n',resp) || strcmpi('no',resp),
           return
        end
        fprintf('Erasing t-test results stored with these data.\n');
    end
    GND.t_tests=[];
end

[n_chan, n_tpt, n_bin, n_sub]=size(GND.indiv_erps);
if strcmpi(method,'boxcar')
    if VERBLEVEL>1,
       fprintf('Downsampling data by a factor of %d after low-pass filtering with a boxcar moving average.\n', ...
           decfactor);
    end
    
    %decimate with moving boxcar window of length decfactor*2-1
    neo_n_tpt=length(decfactor:decfactor:(n_tpt-decfactor+1));
    
    %Recompute time points
    neo_times=zeros(1,neo_n_tpt);
    ct=0;
    for t=decfactor:decfactor:(n_tpt-decfactor+1),
        ct=ct+1;
        if ~rem(ct,10),
            fprintf('Now computing time point %d of %d.\n',ct,neo_n_tpt);
        end
        tstart=t-decfactor+1;
        tend=t+decfactor-1;
        neo_times(ct)=mean(GND.time_pts(tstart:tend));
        for s=1:n_sub,
            for c=1:n_chan,
                for b=1:n_bin,
                    GND.indiv_erps(c,ct,b,s)=mean(GND.indiv_erps(c,tstart:tend,b,s));
                end
                %Cal pulses
                if ~isempty(GND.cals),
                    GND.cals.indiv_cals(c,ct,s)=mean(GND.cals.indiv_cals(c,tstart:tend,s));
                end
            end            
        end
    end
    GND.indiv_erps=GND.indiv_erps(:,1:neo_n_tpt,:,:);
    if ~isempty(GND.cals),
        GND.cals.indiv_cals= GND.cals.indiv_cals(:,1:neo_n_tpt,:);
    end
    GND.time_pts=neo_times;
else
    if VERBLEVEL>1,
       fprintf('Downsampling data by a factor of %d after low-pass filtering with a MATLAB derived FIR filter.\n', ...
           decfactor);
    end
    
    %Recompute time points
    GND.time_pts=round(decimate(GND.time_pts,decfactor,'fir'));
    neo_n_tpt=length(GND.time_pts);
    
    
    %Decimate individual participant ERPs
    for c=1:n_chan,
        fprintf('Now downsampling Channel %d (%s).\n',c,GND.chanlocs(c).labels);
        for b=1:n_bin,
            for s=1:n_sub,
                GND.indiv_erps(c,1:neo_n_tpt,b,s)=decimate(GND.indiv_erps(c,:,b,s),decfactor,'FIR');
            end
        end
    end
    GND.indiv_erps=GND.indiv_erps(:,1:neo_n_tpt,:,:);
        
    if ~isempty(GND.cals),
        %Decimate cal pulse ERPs
        fprintf('Now downsampling cal pulses.\n');
        for s=1:n_sub,
            for c=1:n_chan,
                GND.cals.indiv_cals(c,1:neo_n_tpt,s)=decimate(GND.cals.indiv_cals(c,:,s),decfactor,'FIR');
            end
        end
        GND.cals.indiv_cals=GND.cals.indiv_cals(:,1:neo_n_tpt,:);
    end
end

%Shrink grand fields for overwriting
GND.grands=GND.grands(:,1:neo_n_tpt,:);
GND.grands_stder=GND.grands_stder(:,1:neo_n_tpt,:);
GND.grands_t=GND.grands_t(:,1:neo_n_tpt,:);

%Rebaseline individual ERPs/cal pulses and recompute grands
GND=baselineGND(GND,bsln_wind);

%Recompute sampling rate
GND.srate=GND.srate/decfactor;

%add comand to history
n_hist=length(GND.history);
if isempty(bsln_wind),
    GND.history{n_hist+1}=sprintf('GND=decimate(GND,%d,''%s'',[],%d);', ...
        decfactor,method,VERBLEVEL);
else
    GND.history{n_hist+1}=sprintf('GND=decimate(GND,%d,''%s'',[%d %d],%d);', ...
        decfactor,method,bsln_wind(1),bsln_wind(2), ...
        VERBLEVEL);
end


GND.saved='no';
if strcmpi(save_GND,'yes'),
    GND=save_matmk(GND,'gui');
end
