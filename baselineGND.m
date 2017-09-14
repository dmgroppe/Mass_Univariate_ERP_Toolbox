% baselineGND() - Baseline the ERPs in a Mass Univariate ERP Toolbox GND 
%                 struct variable by removing the mean amplitude within a 
%                 specified time window.
%
% Usage:
%  >> GND=baselineGND(GND,bsln_wind,verblevel);
%
% Required Inputs:
%   GND - A Mass Univariate ERP Toolbox GND structure variable.  To create a 
%         GND variable from Kutaslab ERP files (e.g., *.mas files) use 
%         avgs2GND.m.  To do the same from EEGLAB *.set files use sets2GND.m.
%         See Mass Univariate ERP Toolbox documentation for detailed 
%         information about the format of a GND variable. 
%
% Optional Inputs:
%   bsln_wind  - [vector or NaN] Two element vector specifying the beginning and
%                end (in ms) of the baseline time window (e.g., [-100 -4])
%                or NaN. If NaN, nothing is done and GND variable is returned
%                unchanged. The mean amplitude across all time points within and 
%                including those times will be removed from each ERP. 
%                {default: all time points before 0}
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
% Outputs:
%   GND           - Mass Univariate ERP Toolbox GND structure variable.  
%                   This is the same as the input GND variable, but individual 
%                   ERPs (including cal pulses) have been baselined and grand
%                   averages have been re-computed.  The field GND.bsln_wind 
%                   will be updated to the new baseline window as well.
%
% Notes:
% -GND variable is NOT saved to disk after baselining
% -Function will erase any t-test results stored with GND
% variable since their results may no longer be accurate
%
% Author:
% David Groppe
% Kutaslab, 3/2010

% Changes:
% 8/23/2012 - NaN now a possible value for bsln_wind to avoid any
% baselining. GND is not changed at all when NaN is the baseline window.
%
% 1/10/2017 - ERP grand averages and t-scores are still computed when NaN
% is the baseline window.

function GND=baselineGND(GND,bsln_wind,verblevel)

if nargin<2
    bsln_wind=[];
elseif ~isempty(bsln_wind),
    if isnan(bsln_wind)
        fprintf('Not baselining data.\n');
    else
        if length(bsln_wind)~=2,
            error('Argument bsln_wind needs to be a two element vector.');
        end
        if bsln_wind(2)<bsln_wind(1),
            error('First value of bsln_wind needs to be less than or equal to second value.');
        end
    end
end

global VERBLEVEL;
if nargin<3
    if isempty(VERBLEVEL),
        VERBLEVEL=2;
    end
else
    VERBLEVEL=verblevel;
end

%Erase any t-results as they may not be valid anymore
if ~isempty(GND.t_tests),
    %but ask first
    if VERBLEVEL>1,
        resp=[];
        while ~strcmpi(resp,'y') && ~strcmpi(resp,'n') && ~strcmpi(resp,'yes') ...
                && ~strcmpi(resp,'no'),
            fprintf('These data have t-test results stored with them that may not be accurate after the data have been baselined.\nFor this reason, they will be erased.\n');
            resp=input(sprintf('Continue with baselining (t-test results will be erased)? [y or n] '),'s');
        end
        if strcmpi('n',resp) || strcmpi('no',resp),
           return
        end
        fprintf('Erasing t-test results stored with these data.\n');
    end
    GND.t_tests=[];
end

%find baseline window time points
[n_chan, n_tpt, n_bin, n_sub]=size(GND.indiv_erps);
if isempty(bsln_wind) || sum(isnan(bsln_wind))==0,
    if ~isempty(bsln_wind),
        bsln_tpt(1)=find_tpt(bsln_wind(1),GND.time_pts);
        bsln_tpt(2)=find_tpt(bsln_wind(2),GND.time_pts);
        if VERBLEVEL>1,
            fprintf('Baselining data from %d to %d ms (that''s time point %d to %d).\n', ...
                GND.time_pts(bsln_tpt(1)),GND.time_pts(bsln_tpt(2)),bsln_tpt(1),bsln_tpt(2));
        end
    else
        bsln_tpt=[];
        if VERBLEVEL,
            ids=find(GND.time_pts<0);
            if isempty(ids),
                error('Attempted to use default baseline of all time points before 0. However, all time points in these data are after 0.');
            else
                bsln_tpt(1)=ids(1);
                bsln_tpt(2)=ids(end);
                if VERBLEVEL>1,
                    fprintf('Using default baseline of all time points before 0.\n');
                    fprintf('That''s from %d to %d ms (time points %d to %d).\n', ...
                        GND.time_pts(bsln_tpt(1)),GND.time_pts(bsln_tpt(2)), ...
                        bsln_tpt(1),bsln_tpt(2));
                end
            end
        end
    end
    
    %baseline individual ERPs
    bsln_tpts=bsln_tpt(1):bsln_tpt(2);
    for s=1:n_sub,
        if VERBLEVEL>1,
            fprintf('Baselining data from Participant #%d (%s).\n',s,GND.indiv_subnames{s});
        end
        for b=1:n_bin,
            GND.indiv_erps(:,:,:,s)=reshape(rmbase(GND.indiv_erps(:,:,:,s), ...
                n_tpt,bsln_tpts),n_chan,n_tpt,n_bin);
        end
        if ~isempty(GND.cals),
            cal_size=size(GND.cals.indiv_cals);
            GND.cals.indiv_cals=reshape(rmbase(GND.cals.indiv_cals, ...
                cal_size(2),bsln_tpts),cal_size(1),cal_size(2),n_sub);
        end
    end
end

%Recompute grands, stders, & grand t-scores
for b=1:n_bin,
    bin_subs=find(GND.indiv_bin_ct(:,b));
    GND.sub_ct(b)=length(bin_subs);
    if GND.sub_ct(b),
        GND.grands(:,:,b)=mean(GND.indiv_erps(:,:,b,bin_subs),4);
        GND.grands_stder(:,:,b)=std(GND.indiv_erps(:,:,b,bin_subs),0,4)/sqrt(GND.sub_ct(b));
        GND.grands_t(:,:,b)=GND.grands(:,:,b)./GND.grands_stder(:,:,b);
    else
        watchit(sprintf('No average files contribute to bin %d.',b));
    end
end
%Recompute grand cal pulses
if ~isempty(GND.cals),
    GND.cals.grand_cals=mean(GND.cals.indiv_cals,3);
end

if isnan(bsln_wind),
    GND.bsln_wind=NaN;
else
    GND.bsln_wind=[GND.time_pts(bsln_tpt(1)) GND.time_pts(bsln_tpt(2))];
end

GND.saved='no';


