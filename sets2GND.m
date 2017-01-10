% sets2GND() - Create a Mass Univariate ERP Toolbox GND struct variable from a set of
%              EEGLAB *.set files
%
% Usage:
%  >> GND=sets2GND(gui_infiles_or_tmplt,varargin);
%
% Required Inputs:
%   gui_infiles_or_tmplt - ['gui', a cell array of strings, or a string template]
%                          If 'gui', a GUI is created that allows you to select
%                          which set files to import (this is probably the
%                          easiest way to import files).  If a cell array of
%                          of strings, each element of the cell array should
%                          contain the filename of an EEGLAB set file (e.g.,
%                          {'visodbl01.set','visodbl02.set'}.  Otherwise, this
%                          input should be a filename template (i.e., a string with
%                          # where the subject number should be--'visodbl#.set').
%                          If you provide a template, you must use the option
%                          'sub_ids' (see below) to specify the ID numbers of
%                          each subject. Include the files' path unless the files
%                          are in the current working directory.
%
%
% Optional Inputs:
%   sub_ids          - [integer vector] A set of integers specifying the
%                      subject ID numbers to include in the grand average.
%                      Only necessary if a filename template is given as
%                      the input to gui_infiles_or_tmplt.
%
%   bsln             - [vector or NaN] A pair of numbers (in milliseconds)
%                      specifying the start and end times of the ERP baseline
%                      window or NaN.  If NaN, data are not baselined.  
%                      Otherwise, the mean voltage in the baseline window 
%                      will be removed from each ERP. {default: use all time
%                      points before time 0}.
%
%   use_bins         - [integer vector] A set of integers specifying which
%                      bins to import into MATLAB.  If not specified, all
%                      bins will be imported. Note, if you import only a
%                      subset of bins, the bin numberings will start at 1
%                      and run to the number you've imported (i.e., they may
%                      differ from the bin numbers in the set files.
%                      {default: import all bins}
%
%   exclude_chans    - A cell array of channel labels to exclude from the
%                      importation (e.g., {'A2','lle','rhe'}). You cannot
%                      use both this option and 'include_chans' (below).{default:
%                      not used, import all channels}
%
%   include_chans    - A cell array specifying a subset of channel labels to import
%                      (e.g., {'Fz','Cz','Pz'}).  All other channels will
%                      be ignored. You cannot use both this option and
%                      'exclude_chans' (above). {default: not used, import
%                      all channels}
%
%   exp_name         - [string] Name of the experiment. {default: 'An
%                      Experiment'}
%
%   out_fname        - [string] Filename to save GND variable to.  If empty
%                      (i.e., not specified), a GUI will be created to prompt
%                      you for a filename. If the string 'no save', the GND
%                      variable will not be saved to disk. If no file
%                      extension is given '.GND' will be added to the
%                      filename. {default: not specified}
%
%   verblevel        - An integer specifiying the amount of information you want
%                      this function to provide about what it is doing during runtime.
%                       Options are:
%                        0 - quiet, only show errors, warnings, and EEGLAB reports
%                        1 - stuff anyone should probably know
%                        2 - stuff you should know the first time you start working
%                            with a data set {default value}
%                        3 - stuff that might help you debug (show all
%                            reports)
%
% Output:
%   GND  - Struct variable containing grand averages, individual subject
%          ERPs, mass univariate t-test results and more.
%
%
% Global Variables:
%   VERBLEVEL - Mass Univariate ERP Toolbox level of verbosity (i.e., tells 
%               functions how much to report about what they're doing during
%               runtime) set by the optional function argument 'verblevel'
%
% Notes:
% -This function expects the EEG variable to have information about bins
% stored in it.  More specifically, there needs to be (1) an EEG.bindesc
% field that contains a text description of each bin and (2) epochs of data
% labeled as bin types (e.g., 'bin1').  Bin types are stored in the field
% EEG.epoch(#).eventtype.  Use the function bin_info2EEG.m to add bin
% information to an EEG variable.
%
% -ICs labeled as artifacts will be removed from the data before ERPs are
% formed.  This is done by the function remove_artifact_ics.m.  IC labeling
% can be done via EEGLAB or MATLABmk functions/conventions.  If both
% conventions appear in the EEG variable, only MATLABmk IC labels will be used.
%
% -This function is not yet compatible with ERPLAB's conventions for
% storing bin information.
%
% -The GND fields odelay, cals, indiv_bin_raw_ct, and condesc are specific 
% to Kutaslab data. Other labs should be able to ignore them (though 
% indiv_bin_raw_ct is used by ERPLAB processed data as well.
%
% Author:
% David Groppe
% Kutaslab, 7/2010

%%%%%%%%%%%%%%%% Revision History  %%%%%%%%%%%%%%%%%
% 10/31/2010: 'use_bin' option fixed and problem with reshaping data when no
% ICs are labeled as artifacts fixed-TPU
%
% 11/2/2010: Can now deal with condition codes for Kutaslab data (before it
% assumed that all bins had a condition code of 1)
%
% 11/3/2010: Function now checks to make sure GND file can be written to
% disk before attempting to do so.  Code courtesy of Tom Urbach.
%
% 11/14/2010: 'include_chans' option now reports if function could not find
% a channel the user wished to include. 'include' & 'exclude' options made
% case insensitive
%
% 12/9/2010: 'rawtrials_per_bin' recorded keeping wasn't working if only a
% subset of possible bins were imported.  Works now.
%
% 12/15/2010: 'locfiles' optional input removed.  It wasn't doing anything.
%
% 4/4/2011: Made compatible with Windows.
%
% 8/23/2012: NaN now a possible value for bsln to avoid any baselining
%
% 3/12/2013: If you provide a GND variable filename, program no longer
% asks user if s/he wants to save the file. Thanks to Aaron Newman for the
% suggestion.

%%%%%%%%%%%%%%%% Future Work  %%%%%%%%%%%%%%%%%
%
% Could make function able to read ERPLAB bin information in set files 
%
% Add odelay for Kutaslab?  Currently GND.odelay=[];
%
% GND.indiv_traits=[]; <-make import option for this (i.e., function
% should be able to import this info from text file)
%
% Make import option for GND.indiv_traits this (i.e., sets2GND.m should be able to import
% this info from text file)

function GND=sets2GND(gui_infiles_or_tmplt,varargin)

p=inputParser;
p.addRequired('gui_infiles_or_tmplt',@(x) ischar(x) || iscell(x));
p.addParamValue('sub_ids',[],@isnumeric);
p.addParamValue('bsln',[],@(x) sum(isnan(x)) || (isnumeric(x) && length(x)==2));
p.addParamValue('exp_name','An Experiment',@ischar);
p.addParamValue('use_bins',[],@isnumeric);
p.addParamValue('exclude_chans',[],@(x) ischar(x) || iscell(x));
p.addParamValue('include_chans',[],@(x) ischar(x) || iscell(x));
p.addParamValue('out_fname',[],@ischar);
p.addParamValue('verblevel',[],@(x) isnumeric(x) && (length(x)==1));
p.parse(gui_infiles_or_tmplt,varargin{:});


global EEG
global VERBLEVEL;

if isempty(p.Results.verblevel),
    if isempty(VERBLEVEL),
        VERBLEVEL=2;
    end
else
    VERBLEVEL=p.Results.verblevel;
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


%%  Select files for loading
if strcmpi(gui_infiles_or_tmplt,'GUI'),
    loading=1;
    infiles=[];
    while loading,
        [neofname, inpath]=uigetfile({'*.set','*.set files'; ...
            '*.*','All Files (*.*)'},'EEGLAB set Files to Import','MultiSelect','on');
        if ischar(neofname),
            clear infname;
            infname{1}=neofname; %make it a cell array for consistent syntax below
        else
            infname=neofname;
        end
        if ~inpath,
            if isempty(infiles),
                fprintf('File selection cancelled.  Aborting sets2GND.\n');
                GND=[];
                return;
            else
                loading=0;
            end
        else
            if isempty(infiles),
                infiles=cell(1,length(infname)); %preallocate mem
                for a=1:length(infname),
                    infiles{a}=[inpath infname{a}];
                end
            else
                n_files=length(infiles);
                ct=0;
                for a=(n_files+1):(n_files+length(infname)),
                    ct=ct+1;
                    infiles{a}=[inpath infname{ct}];
                end
            end
            % trigger GUI to see if user wants to load more
            resp=questdlg('Do you want to load more files?',...
                'OUTPUT','Yes','No','Yes');
            if strcmpi(resp,'No'),
                loading=0;
            end
        end
    end
elseif iscell(p.Results.gui_infiles_or_tmplt)
    infiles=p.Results.gui_infiles_or_tmplt;
else
    if isempty(p.Results.sub_ids)
        error('If gui_infiles_or_tmplt is a string filename template, you must specify the subject numbers with the argument ''sub_ids''.');
    elseif ~ismember('#',p.Results.gui_infiles_or_tmplt),
        error('The filename template %s needs to contain a # to indicate where the subject numbers go (e.g., odbl#.nrm).', ...
            p.Results.gui_infiles_or_tmplt);
    else
        lb_id=find(p.Results.gui_infiles_or_tmplt=='#');
        prefix=p.Results.gui_infiles_or_tmplt(1:(lb_id(1)-1)); %indexing lb_id by 1 in case there are multiple pound signs
        postfix=p.Results.gui_infiles_or_tmplt((lb_id(1)+1):end);
        n_infiles=length(p.Results.sub_ids);
        infiles=cell(1,n_infiles);
        for s=1:n_infiles,
            no_pad=[prefix num2str(p.Results.sub_ids(s)) postfix];
            if p.Results.sub_ids(s)<10,
                padded=[prefix num2str(p.Results.sub_ids(s),'%.2d') postfix];
                if ismac || isunix
                    %Mac or Unix OS
                    [sP wP]=unix(['ls ' padded]); %if sp==0, then the file exists, need wP argument so that it isn't automatically displayed in command window
                    [sNP wNP]=unix(['ls ' no_pad]); %Need wP argument so that it isn't automatically displayed in command window
                    if ~sP
                        if ~sNP,
                            error('You have a file named %s and one named %s.  I don''t know which one to load!\n', ...
                                padded,no_pad);
                        else
                            infiles{s}=padded;
                        end
                    else
                        infiles{s}=no_pad;
                    end
                else
                    %Must be a PC
                    infiles{s}=no_pad;
                end
            else
                infiles{s}=no_pad;
            end
        end
    end
end

%Get rid of any redundant input files
infiles=unique(infiles);
n_infiles=length(infiles);


%% Load first data set and draft GND fields
EEG=pop_loadset(infiles{1});
EEG_var_check(EEG,infiles{1});

[n_EEG_chans, n_pts, n_epochs]=size(EEG.data);
if ~isempty(exclude_chans),
    use_chans=zeros(1,n_EEG_chans);
    ex_ct=0;
    for a=1:n_EEG_chans,
        if ~ismember_ci(EEG.chanlocs(a).labels,exclude_chans)
            use_chans(a)=1;
        else
            ex_ct=ex_ct+1;
            ex_labels{ex_ct}=EEG.chanlocs(a).labels;
        end
    end
    n_chans=sum(use_chans);
    use_chans=find(use_chans==1);
    
    if VERBLEVEL>1,
        missed=setdiff_ci(exclude_chans,ex_labels);
        n_missed=length(missed);
        if n_missed,
            if n_missed==1,
                msg=sprintf('I attempted to exclude the following channel, but it was not found:');
            else
                msg=sprintf('I attempted to exclude the following channels, but they were not found:');
            end
            for a=1:n_missed,
                msg=[msg ' ' missed{a}];
            end
            watchit(msg);
        end
    end
elseif ~isempty(include_chans),
    use_chans=zeros(1,n_EEG_chans);
    in_ct=0;
    for a=1:n_EEG_chans,
        if ismember_ci(EEG.chanlocs(a).labels,include_chans)
            use_chans(a)=1;
            in_ct=in_ct+1;
            in_labels{in_ct}=EEG.chanlocs(a).labels;
        end
    end
    n_chans=sum(use_chans);
    use_chans=find(use_chans==1);
    
    if VERBLEVEL>1,
        missed=setdiff_ci(include_chans,in_labels);
        n_missed=length(missed);
        if n_missed,
            if n_missed==1,
                msg=sprintf('I attempted to include the following channel, but it was not found:');
            else
                msg=sprintf('I attempted to include the following channels, but they were not found:');
            end
            for a=1:n_missed,
                msg=[msg ' ' missed{a}];
            end
            watchit(msg);
        end
    end
else
    n_chans=n_EEG_chans;
    use_chans=1:n_chans;
end

n_EEG_bins=length(EEG.bindesc);
if isempty(p.Results.use_bins),
    use_bins=1:n_EEG_bins;
    n_bins=n_EEG_bins;
else
    %modicum of error checking
    if min(p.Results.use_bins)<1,
        error('You cannot request to import bins with indices less than 1 (e.g., "Bin -1").');
    elseif length(unique(p.Results.use_bins)) ~= length(p.Results.use_bins)  % TPU checking for duplicate bin numbers in use_bins
        error('You cannot request to import duplicate bins\n');              % TPU
    elseif max(p.Results.use_bins)>n_EEG_bins,
        error('File %s only has %d bins, but you''ve requested to import Bin %d.\n', ...
            infiles{1},n_EEG_bins,max(p.Results.use_bins));
    end
    use_bins=p.Results.use_bins;
    n_bins=length(use_bins);
end
%initialize cell array for keeping track of used urevents
used_urevents=cell(n_infiles,n_bins); %at most there will be n_infiles subjects

GND.exp_desc=p.Results.exp_name;
GND.filename=[];
GND.filepath=[];
GND.saved='no';
GND.grands=zeros(n_chans,n_pts,n_bins)*NaN;
GND.grands_stder=GND.grands;
GND.grands_t=GND.grands;
GND.sub_ct=zeros(1,n_bins);
GND.chanlocs=EEG.chanlocs(use_chans);
GND.bin_info=[];
if isfield(EEG,'condesc'),
    for b=1:n_bins,
        GND.bin_info(b).bindesc=EEG.bindesc{use_bins(b)};
        GND.bin_info(b).condcode=EEG.binccodes(use_bins(b));
    end
    GND.condesc=EEG.condesc;
else
    for b=1:n_bins,
        GND.bin_info(b).bindesc=EEG.bindesc{use_bins(b)};
        GND.bin_info(b).condcode=1; %Not Kutaslab data, so condition codes aren't important
    end
    
    %Note, Kutaslab data won't have condition descriptors
    GND.condesc{1}='Experiment (not cal pulses)';
end
GND.time_pts=EEG.times;
if isempty(p.Results.bsln)
    %use all time points before 0 or first time point
    GND.bsln_wind(1)=EEG.times(1);
    ids=find(EEG.times<0);
    if isempty(ids)
        GND.bsln_wind(2)=EEG.times(1);
    else
        GND.bsln_wind(2)=EEG.times(max(ids));
    end
else
    GND.bsln_wind=p.Results.bsln;
end
GND.odelay=[];
GND.srate=EEG.srate;
GND.indiv_fnames=infiles;
if isempty(EEG.subject),
    if VERBLEVEL,
        fprintf('Set file, %s, does not have a subject name. I will use the filename as the subject''s name.\n', ...
            infiles{1});
    end
    GND.indiv_subnames{1}=infiles{1};
else
    GND.indiv_subnames{1}=EEG.subject;
end
GND.indiv_traits=[];
GND.indiv_bin_ct=zeros(n_infiles,n_bins);
if isfield(EEG,'rawtrials_per_bin'),
    GND.indiv_bin_raw_ct=zeros(n_infiles,n_bins);
else
    GND.indiv_bin_raw_ct=zeros(n_infiles,n_bins)*NaN; %This field is only relevant to Kutaslab data.  It keeps track of how many trials there were before artifact rejection.
end
GND.indiv_erps=zeros(n_chans,n_pts,n_bins,n_infiles);
GND.indiv_art_ics=[];
if isfield(EEG,'cal_info'),
    %This field is only relevant to Kutaslab data.
    cal_size=size(EEG.cal_info.erps);
    GND.cals.indiv_cals=zeros(n_chans,cal_size(2),n_infiles);
    GND.cals.indiv_cal_ct=zeros(1,n_infiles);
    GND.cals.grand_cals=zeros(cal_size);
    GND.cals.caldesc='cal pulses';
    GND.cals.condcode=0;
    GND.cals.condesc='calibration';
    cal_info_present=1;
else
    cal_info_present=0;
    GND.cals=[];
end
GND.history=[];
GND.t_tests=[];

%% Loop through set files
sub_ct=1;
for filenum=1:n_infiles,
    new_sub=1; %Assume data from this subject has not been already loaded, until we learn otherwise (This is to deal with the fact that data from the same participant may be distributed among multiple set files)
    if filenum>1,
        global EEG; %this needs to be done to be compatible with EEGLAB (which stores EEG globally)
        if VERBLEVEL,
            fprintf('\n\n');
        end
        EEG=pop_loadset(infiles{filenum});
        
        %collect sub names
        if isempty(EEG.subject),
            if VERBLEVEL,
                fprintf('Set file, %s, does not have a subject name. I will use the filename as the subject''s name.\n', ...
                    infiles{filenum});
            end
            sub_name=infiles{filenum};
        else
            sub_name=EEG.subject;
        end
        ur_sub_id=0;
        for old_sub=1:sub_ct,
            if strcmpi(GND.indiv_subnames{old_sub},sub_name),
                ur_sub_id=old_sub;
                if VERBLEVEL
                    fprintf('Set file, %s, will be appended to data already loaded from Subject %s.\n', ...
                        infiles{filenum},GND.indiv_subnames{old_sub});
                end
                new_sub=0;
                break;
            end
        end
        if ~ur_sub_id,
            sub_ct=sub_ct+1;
            sub=sub_ct;
            GND.indiv_subnames{sub}=sub_name;
        else
            sub=ur_sub_id;
        end
        
        EEG_var_check(EEG,infiles{filenum});
        
        [n_EEG_chans2, n_pts2, n_epochs]=size(EEG.data);
        n_EEG_bins2=length(EEG.bindesc);
        
        %Figure out which channels to use
        if ~isempty(exclude_chans),
            use_chans=zeros(1,n_EEG_chans2);
            ex_ct=0;
            ex_labels=[];
            for a=1:n_EEG_chans2,
                if ~ismember_ci(EEG.chanlocs(a).labels,exclude_chans)
                    use_chans(a)=1;
                else
                    ex_ct=ex_ct+1;
                    ex_labels{ex_ct}=EEG.chanlocs(a).labels;
                end
            end
            n_chans2=sum(use_chans);
            use_chans=find(use_chans==1);
            
            if VERBLEVEL>1,
                missed=setdiff_ci(exclude_chans,ex_labels);
                n_missed=length(missed);
                if n_missed,
                    if n_missed==1,
                        msg=sprintf('I attempted to exclude the following channel, but it was not found:');
                    else
                        msg=sprintf('I attempted to exclude the following channels, but they were not found:');
                    end
                    for a=1:n_missed,
                        msg=[msg ' ' missed{a}];
                    end
                    watchit(msg);
                end
            end
        elseif ~isempty(include_chans),
            use_chans=zeros(1,n_EEG_chans2);
            in_ct=0;
            in_labels=[];
            for a=1:n_EEG_chans2,
                if ismember_ci(EEG.chanlocs(a).labels,include_chans)
                    use_chans(a)=1;
                    in_ct=in_ct+1;
                    in_labels{in_ct}=EEG.chanlocs(a).labels;
                end
            end
            n_chans2=sum(use_chans);
            use_chans=find(use_chans==1);
            
            if VERBLEVEL>1,
                missed=setdiff_ci(include_chans,in_labels);
                n_missed=length(missed);
                if n_missed,
                    if n_missed==1,
                        msg=sprintf('I attempted to include the following channel, but it was not found:');
                    else
                        msg=sprintf('I attempted to include the following channels, but they were not found:');
                    end
                    for a=1:n_missed,
                        msg=[msg ' ' missed{a}];
                    end
                    watchit(msg);
                end
            end
        else
            n_chans2=n_EEG_chans2;
            use_chans=1:n_chans2;
        end
        
        %Check for consistency with GND
        if EEG.srate~=GND.srate,
            error('File %s has a different sampling rate than file %s.\n', ...
                infiles{1},infiles{filenum});
        end
        if n_chans2~=n_chans,
            error('File %s has a different number of channels to import than file %s.\n', ...
                infiles{1},infiles{filenum});
        end
        if n_pts2~=n_pts,
            error('File %s has a different number of time points than file %s.\n', ...
                infiles{1},infiles{filenum});
        end
        if ~min(EEG.times==GND.time_pts),
            error('The epochs in file %s begin and end at different times than file %s.\n', ...
                infiles{1},infiles{filenum});
        end
        if isempty(p.Results.use_bins),
            if n_EEG_bins2~=n_bins,
                error('File %s has a different number of bins than file %s.\n', ...
                    infiles{1},infiles{filenum});
            end
        else
            if max(p.Results.use_bins)>n_EEG_bins2,
                error('File %s only has %d bins, but you''ve requested to import Bin %d.\n', ...
                    infiles{filenum},n_EEG_bins2,max(p.Results.use_bins));
            end
        end
        if isfield(EEG,'condesc'),
            newfile_condesc=EEG.condesc;
        else
            newfile_condesc{1}='Experiment (not cal pulses)';
        end
        if ~isequal(GND.condesc,newfile_condesc),
            error('The condition code descriptor (EEG.condesc) in file %s is not the same as that in previous files.\n', ...
                infiles{filenum});
        end
        bin_ct=0;
        for b=use_bins,
            bin_ct=bin_ct+1;
            if ~strcmpi(GND.bin_info(bin_ct).bindesc,EEG.bindesc{b}),
                error('The #%d imported bin in file %s is different than that of file %s.\n', ...
                    b,infiles{1},infiles{filenum});
            end
        end
        try
            [fs1, fs2, er]=comp_struct_quiet(GND.chanlocs,EEG.chanlocs(use_chans));
        catch
            error('File %s''s imported channel location information differs from that of file %s.\n', ...
                infiles{1},infiles{filenum});
        end
        if ~isempty(er),
            error('File %s''s imported channel location information differs from that of file %s.\n', ...
                infiles{1},infiles{filenum});
        end
    else
        sub=sub_ct; %for first set file, sub=1 (since sub_ct=1)
    end
    
    %Remove ICs labeled as artifacts
    fldnames=fieldnames(EEG);
    art_ics=[];
    ics_removed=0;
    if sum(EEG.reject.gcompreject),
        ics_removed=1;
    else
        for fn=1:length(fldnames),
            if strcmpi(fldnames{fn},'iclabels'),
                ics_removed=1;
                break;
            end
        end
    end
    if ics_removed,
        art_ics=remove_artifact_ics(GND.bsln_wind);
    end
    GND.indiv_art_ics{filenum}=art_ics;
    
    %baseline data (if not already baselined by ics_removed)
    if isnan(GND.bsln_wind)
        fprintf('NOT baselining data.\n');
    elseif ~isempty(GND.bsln_wind) && ~ics_removed,
        bsln_pts(1)=find_tpt(GND.bsln_wind(1),EEG.times);
        bsln_pts(2)=find_tpt(GND.bsln_wind(2),EEG.times);
        if VERBLEVEL>=2,
            fprintf('Baselining data by removing mean EEG between %d and %d ms (time points %d and %d).\n', ...
                GND.bsln_wind(1),GND.bsln_wind(2),bsln_pts(1),bsln_pts(2));
        end
        
        EEG.data=rmbase(EEG.data,n_pts,bsln_pts(1):bsln_pts(2));
        EEG.data=reshape(EEG.data,n_EEG_chans,n_pts,n_epochs);
    end
    
    
    %Compute ERPs, use GND.indiv_erps as a running sum
    for a=1:n_epochs,
        for b=1:length(EEG.epoch(a).eventtype),
            if (length(EEG.epoch(a).eventtype{b})>2) && strcmpi(EEG.epoch(a).eventtype{b}(1:3),'bin')
                bin=str2num(EEG.epoch(a).eventtype{b}(4:end));
                if ismember(bin,use_bins) && ~index_cellarray_or_vector(EEG.epoch(a).eventlatency,b,VERBLEVEL),
                    
                    % BEGIN TPU
                    binind=find(use_bins==bin);  % lookup dimension of array to load
                    if (length(binind)>1) % die if duplicate bins ... this should be caught earlier
                        error('ERROR sets2GND.m Duplicate bins in use_bins');
                    end
                    % END TPU
                    
                    %only consider events with epoch latency of 0
                    new_event=1;
                    urevent=index_cellarray_or_vector(EEG.epoch(a).eventurevent,b,VERBLEVEL);
                    if ~new_sub,
                        %check to see if event has already been used for
                        %this bin
                        if ismember(urevent,used_urevents{sub,binind}), % TPU
                            if VERBLEVEL>1,
                                fprintf('Epoch %d has already been added to Bin %d for Subject %s.\n', ...
                                    a,bin,GND.indiv_subnames{sub});
                            end
                            new_event=0;
                        end
                    end
                    if new_event,
                        GND.indiv_erps(:,:,binind,sub)=GND.indiv_erps(:,:,binind,sub)+EEG.data(use_chans,:,a); % TPU bin->binind
                        GND.indiv_bin_ct(sub,binind)=GND.indiv_bin_ct(sub,binind)+1;   % TPU bin->binind
                        used_urevents{sub,binind}=[used_urevents{sub,binind} urevent]; % TPU bin->binind
                    end
                end
            end
        end
    end
    
    if isfield(EEG,'rawtrials_per_bin'),
        %Again, this field is only relevant to Kutaslab data.  It keeps
        %track of how many trials there were before artifact rejection.
        GND.indiv_bin_raw_ct(sub,:)=GND.indiv_bin_raw_ct(sub,:)+EEG.rawtrials_per_bin(use_bins);
    end
    
    if isfield(EEG,'cal_info')
        %Again, this field is only relevant to Kutaslab data.
        %Note, I assume that even if a single participant's data are
        %distributed among multiple set files, EEG.cal_info is the same in
        %all the set files.
        GND.cals.indiv_cals(use_chans,:,sub)=EEG.cal_info.erps(use_chans,:);
        GND.cals.indiv_cal_ct(sub)=EEG.cal_info.npulse_used;
        GND.cals.caldesc='cal pulses';
        GND.cals.condcode=0;
        GND.cals.condesc='calibration';
        cal_info_present=1;
    end
    
    clear global EEG; %make sure it's not saved since ICs might have been removed
end %filenum loop


%remove unused elements of GND variable (if any), because there are fewer
%subjects than set files
GND.indiv_bin_ct=GND.indiv_bin_ct(1:sub_ct,:);
GND.indiv_bin_raw_ct=GND.indiv_bin_raw_ct(1:sub_ct,:);
GND.indiv_erps=GND.indiv_erps(:,:,:,1:sub_ct);
if cal_info_present,
    GND.cals.indiv_cals(:,:,1:sub_ct)=GND.cals.indiv_cals(:,:,1:sub_ct);
    GND.cals.indiv_cal_ct=GND.cals.indiv_cal_ct(1:sub_ct);
end

%turn running sums into means
for sub=1:sub_ct,
    if VERBLEVEL>1
        fprintf('\nTrials per bin for Subject %s:\n',GND.indiv_subnames{sub});
    end
    for b=1:n_bins,
        bin_ct=GND.indiv_bin_ct(sub,b);
        if bin_ct,
            if VERBLEVEL>1,
                fprintf('Bin %d (%s): %d trials\n',b,GND.bin_info(b).bindesc,bin_ct);
            end
            GND.indiv_erps(:,:,b,sub)=GND.indiv_erps(:,:,b,sub)/bin_ct;
            GND.sub_ct(b)=GND.sub_ct(b)+1;
        else
            watchit(sprintf('Subject %s has no epochs that fall into bin %d.',GND.indiv_subnames{sub},b));
            GND.indiv_erps(:,:,b,sub)=GND.indiv_erps(:,:,b,sub)*NaN;
        end
    end
end

%Baseline individual ERPs/cal pulses and compute grands
if VERBLEVEL>1,
    fprintf('\n\n'); %add a couple lines between between counts and baselining info
end
GND=baselineGND(GND,GND.bsln_wind);% This line actually computes the grand average ERPs. It doesn't just baseline them.


%Compute grand average cal pulses (if individual subject cal pulses present)
if cal_info_present,
    GND.cals.grand_cals=mean(GND.cals.indiv_cals,3);
end

%% Save GND variable
if isempty(p.Results.out_fname),
    %Create GUI
    [jname, jpath]=uiputfile({'*.GND','*.GND files'; ...
        '*','All files'},'Save GND variable as:','untitled.GND');
    if ~jpath,
        fprintf('Output filename selection cancelled.  GND variable NOT saved to disk.\n');
    else
        %test to make sure file can be created
        isW=isWriteable([jpath jname]);
        if isW,
            GND=save_matmk(GND,jname,jpath,1); % 1 means that user won't be asked again about saving file
        else
            fprintf('GND file could not be saved, but should still exist in MATLAB workspace.');
        end
    end
elseif ~strcmpi(p.Results.out_fname,'no save'),
    [jpath, jname]=pathNname(p.Results.out_fname);
    %Add .GND extension if no extension given
    if ~ismember('.',jname),
        jname=[jname '.GND'];
    end
    fprintf('Saving GND variable as %s\n',[jpath jname]);
    force=1;
    GND=save_matmk(GND,jname,jpath,force);
end

function EEG_var_check(EEG,fname)

if isempty(EEG.epoch),
    error('File %s stores continuous EEG data. You need to convert the data into epochs.',fname);
elseif ~isfield(EEG,'bindesc'),
    error('The variable stored in file %s does not contain a bin descriptor. You need to apply bin_info2EEG.m to it.',fname);
end

function value=index_cellarray_or_vector(cellarray_or_vector,index,VERBLEVEL)

if index>length(cellarray_or_vector),
    %This if statement is necessary for compatibility with Kutaslab .set
    %files, which can have multiple types per epoch but typically only a
    %single latency value
    index=length(cellarray_or_vector);
    if VERBLEVEL>2,
        fprintf('EEG.epoch(#).eventtype index exceeds EEG.epoch(#).latency index.  I will simply use the last element of EEG.epoch(#).latency.\n');
    end
end

if iscell(cellarray_or_vector),
    value=cellarray_or_vector{index};
else
    value=cellarray_or_vector(index);
end

%%%%%%%% END OF MAIN FUNCTION %%%%%%%


function isW = isWriteable(inFile)
% Function checks to make sure "inFile" can be written to disk. "inFile" is
% not modified by the function.
%
% Author: Tom Urbach
%

isW = 0;
[fid, message] = fopen(inFile,'a+');
if (fid ~= -1)
    isW = 1;
    fclose(fid);
else
    fprintf('Error opening %s: %s\n', inFile, message);
end



function yesno=ismember_ci(str,str_array)
% function yesno=ismember_ci(str,str_array)
% A case insensitive version of ismember.m but just for strings
%
% Inputs:
%  str       - a string
%  str_array - a cell array of strings
%
% Outputs:
%  yesno - 1 if str is a member of str_array.  0 otherwise.  Comparison is
%          not case sensitive.

yesno=0;
n_str=length(str_array);

for a=1:n_str,
    if strcmpi(str,str_array{a})
        yesno=1;
        break;
    end
end



function dif_str=setdiff_ci(superset,subset)
%function dif_str=setdiff_ci(superset,subset)
%
% Inputs:
%   superset - an cell array of strings
%   subset   - an cell array of strings
%
% Outputs:
%   dif_str - a cell array of the strings that are in superset but NOT
%             subset. Comparison is not case sensitive.
%

n_super=length(superset);
n_sub=length(subset);
dif_ct=0;
dif_str=[];
for a=1:n_super,
    found=0;
    for b=1:n_sub,
        if strcmpi(superset{a},subset{b}),
            found=1;
            break
        end
    end
    if ~found,
        dif_ct=dif_ct+1;
        dif_str{dif_ct}=superset{a};
    end
end