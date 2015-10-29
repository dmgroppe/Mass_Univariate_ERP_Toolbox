% bin_info2EEG() - Adds "bin" information to an EEGLAB EEG variable based on
%                  a "bin list file," a text file that indicates which types 
%                  of events belong to which bins.  See the help comments
%                  for the function read_bin_list_file.m for information on
%                  how to create a bin list file. Note, the EEG variable
%                  can contain epoched OR continuous data.
%
% Usage:
%  >>  EEG=bin_info2EEG(EEG_or_fname,blf_fname,save_fname,verblevel);
%
%
% Required Inputs:
%   EEG_or_fname - EEGLAB EEG struct variable containing epochs of EEG data
%                  (i.e., not continuous EEG) or the filename of a set file
%                  that contains such an EEG struct variable.
%
%   blf_fname    - [string] The filename of the bin list file. Be sure to
%                  include the file's path if the file is not in the current 
%                  working directory.
%
% Optional Inputs:
%   save_fname     - The filename with which you want to save the EEG
%                    variable after the bin information has been added.  If
%                    not specified, the EEG variable will not be saved to
%                    disk (it will simply be returned as the output of the
%                    function). {default: not specified}
%   verblevel      - An integer specifiying the amount of information you want
%                    this function to provide about what it is doing during runtime.
%                     Options are:
%                      0 - quiet, only show errors, warnings, and EEGLAB reports
%                      1 - stuff anyone should probably know
%                      2 - stuff you should know the first time you start working
%                          with a data set {default value if not globally specified}
%                      3 - stuff that might help you debug (show all
%                          reports)
%
% Output:
%   EEG        - EEGLAB EEG struct variable with bin information added.
%                More specifically, there will be (1) an EEG.bindesc
%                field that contains a text description of each bin, (2)
%                events of bin types (e.g., 'bin1') will be added to the field
%                EEG.event, and (3) if the data are epoched, bin types will
%                be stored in EEG.epoch(#).eventtype.
%
% Global Variable:
%   VERBLEVEL - Mass Univariate ERP Toolbox level of verbosity (i.e., tells
%               functions how much to report about what they're doing during
%               runtime) set by the optional function argument 'verblevel'
%
% Notes:
% -If bin information has already been added to an EEG variable, this
% function cannot add to it or overwrite it.  You need to remove it and
% then call this function.
%
% -This function is not yet compatible with ERPLAB's conventions for
% storing bin information.
%
% -EEG.urevent field is not modified by this function.  Thus "bin" events
% only are recorded in EEG.event and EEG.epoch fields.
%
% Author:
% David Groppe
% Kutaslab, 9/2010

%%%%%%%%%%%%%%%% Revision History  %%%%%%%%%%%%%%%%%
% 12/10/2010-Now works on continuous data.
%
% 4/13/2011-Now works on EEG variables that don't have
% EEG.event(x).duration fields
%
% 9/5/2011-New 'bin*' elements of EEG.events had the wrong bin #. Also, the 
% code wouldn't work if EEG.epoch(#).eventtype was numeric.  These
% problems have been fixed. Thanks to Haibo "Seapsy" Zhou for finding these
% problems and their solutions.
%
% 6/25/2012-Function crashed or could copy over incorrect EEG.epoch info
% (e.g. EEG.epoch(x).latency) when run on epoched data with more than one 
% event per epoch. Thanks to Clemens Maidhof for alerting me to this.
%


%%%%%%%%%%%%%%%% Future Work  %%%%%%%%%%%%%%%%%
%-Make compatible with ERPLAB

function EEG=bin_info2EEG(EEG_or_fname,blf_fname,save_fname,verblevel)

global VERBLEVEL;

if nargin<3
    save_fname=[];
end

if nargin<4,
    if isempty(VERBLEVEL),
        VERBLEVEL=2;
    end
else
    VERBLEVEL=verblevel;
end

if ischar(EEG_or_fname),
    EEG=pop_loadset(EEG_or_fname);
else
    EEG=EEG_or_fname;
    clear EEG_or_fname;
end

%Are data epoched or continuous?
if length(size(EEG.data))==2
    continuous=1;
else
    continuous=0;
    n_ep=length(EEG.epoch);
end
n_ev=length(EEG.event);

%Check for bin_information. If found, warn that it needs to be erased
if isfield(EEG,'bindesc'),
    error('This EEG variable already has bin information added to it (i.e., there is a ''bindesc'' field). You need to remove it before using this function.');
else
    if continuous,
        for a=1:n_ev,
            if strcmpi(EEG.event(a).type(1:3),'bin'),
                error('This EEG variable already has bin information added to it (i.e., Event #%d is of type ''%s''). You need to remove it before using this function.', ...
                    a,EEG.event(a).type);
            end
        end
    else
        for a=1:n_ep,
            if iscell(EEG.epoch(a).eventtype),
                n_type=length(EEG.epoch(a).eventtype);
            else
                %if only one event per epoch, a string is used instead
                %of a cell array
                n_type=1;
            end
            
            for b=1:n_type,
                if n_type>1
                    if (length(EEG.epoch(a).eventtype{b})>3) && strcmpi(EEG.epoch(a).eventtype{b}(1:3),'bin'),
                        error('This EEG variable already has bin information added to it (i.e., Epoch #%d has a type of ''%s''). You need to remove it before using this function.', ...
                            a,EEG.epoch(a).eventtype{b});
                    end
                else
                    if (length(EEG.epoch(a).eventtype)>3) && strcmpi(EEG.epoch(a).eventtype(1:3),'bin'),
                        error('This EEG variable already has bin information added to it (i.e., Epoch #%d has a type of ''%s''). You need to remove it before using this function.', ...
                            a,EEG.epoch(a).eventtype);
                    end
                end
            end
        end
    end
end

%% Read bin information from bin list file (a text file)
[bindesc, bintype]=read_bin_list_file(blf_fname,VERBLEVEL);
EEG.bindesc=bindesc;
n_bin=length(bindesc);

bin_ct=zeros(1,n_bin);

%% Add bin information
if continuous,
    ev_ct=n_ev;
    for a=1:n_ev,
        for bin=1:n_bin,
            if ismember(EEG.event(a).type,bintype{bin})
                ev_ct=ev_ct+1;
                EEG.event(ev_ct).latency=EEG.event(a).latency;
                if isfield(EEG.event(a),'duration'),
                    %EEGLAB apparently doesn't always add a duration field
                    EEG.event(ev_ct).duration=EEG.event(a).duration;
                end
                EEG.event(ev_ct).type=['bin' int2str(bin)];
                EEG.event(ev_ct).urevent=EEG.event(a).urevent;
                
                %update bin counter
                bin_ct(bin)=bin_ct(bin)+1;
            end
        end
    end
else
    %epoched data
    
    % If the EEG.epoch field contains numeric values, change them to
    % strings since the rest of the code expects this field to contain only
    % strings.
    num2str_flag=0;
    for a=1:n_ep
        if isnumeric(EEG.epoch(a).eventtype)
            temp=EEG.epoch(a).eventtype;
            EEG.epoch(a).eventtype=[];
            for b=1:length(temp)
                num2str_flag=1;
                EEG.epoch(a).eventtype{b}=num2str(temp(b));
            end
        end
    end
    if num2str_flag,
       if VERBLEVEL>=1,
          fprintf('EEG.epoch(#).eventtype field contained numeric values. These have been changed to strings.'); 
       end
    end
    
    add_ev_ct=0;
    for a=1:n_ep,
        if iscell(EEG.epoch(a).eventtype),
            add_ct=0;
            n_ep_type=length(EEG.epoch(a).eventtype);
            for bin=1:n_bin,
                if ismember(EEG.epoch(a).eventtype,bintype{bin})
                    %convert non cell array fields of EEG.epoch(1) into
                    %cell arrays
                    
                    % Note: EEG.epoch(a).event stays a vector even when
                    % multiple events are present in an epoch
                    
                    %Some versions of EEGLAB apparently don't have the
                    %eventduration field or don't create it automatically
                    if isfield(EEG.epoch(a),'eventduration')
                        evdurflag=1;
                        if ~iscell(EEG.epoch(a).eventduration)
                            temp=EEG.epoch(a).eventduration; %EEG.eventduration
                            EEG.epoch(a).eventduration=cell(1,1);
                            EEG.epoch(a).eventduration{1}=temp;
                        end
                    else
                        evdurflag=0;
                    end
                    
                    if ~iscell(EEG.epoch(a).eventlatency)
                        temp=EEG.epoch(a).eventlatency; %EEG.eventlatency
                        EEG.epoch(a).eventlatency=cell(1,1);
                        EEG.epoch(a).eventlatency{1}=temp;
                    end
                    
                    if ~iscell(EEG.epoch(a).eventtype)
                        temp=EEG.epoch(a).eventtype; %EEG.eventtype
                        EEG.epoch(a).eventtype=cell(1,1);
                        EEG.epoch(a).eventtype{1}=temp;
                    end
                    
                    if ~iscell(EEG.epoch(a).eventurevent)
                        temp=EEG.epoch(a).eventurevent; %EEG.eventurevent
                        EEG.epoch(a).eventurevent=cell(1,1);
                        EEG.epoch(a).eventurevent{1}=temp;
                    end
                    
                    % I am not sure why the block of code below was ever
                    % added.  It appears redundant with block of code below
                    % it.  Am just commenting it out for time being.
                    
                    %Add bin info to EEG.epoch
%                     add_ct=add_ct+1;
%                     add_ev_ct=add_ev_ct+1;
%                     EEG.epoch(a).eventtype{add_ct+n_ep_type}=['bin' num2str(bin)];
%                     %Duplicate other information about event
%                     if evdurflag,
%                         EEG.epoch(a).eventduration{add_ct+n_ep_type}=EEG.epoch(a).eventduration{b};
%                     end
%                     EEG.epoch(a).eventlatency{add_ct+n_ep_type}=EEG.epoch(a).eventlatency{b};
%                     EEG.epoch(a).event(add_ct+n_ep_type)=n_ev+add_ev_ct;
%                     EEG.epoch(a).eventurevent{add_ct+n_ep_type}=EEG.epoch(a).eventurevent{b};
%                     
%                     source_ev_id=EEG.epoch(a).event(b);
%                     %Note, it looks like EEGLAB adds an event for a urevent everytime it occurs in an epoch
%                     
%                     %Add bin info to EEG.event
%                     EEG.event(n_ev+add_ev_ct).type=['bin' num2str(bin)];
%                     EEG.event(n_ev+add_ev_ct).latency=EEG.event(source_ev_id).latency;
%                     if isfield(EEG.event(source_ev_id),'duration')
%                         %EEGLAB apparently doesn't always create this field
%                         EEG.event(n_ev+add_ev_ct).duration=EEG.event(source_ev_id).duration;
%                     end
%                     EEG.event(n_ev+add_ev_ct).urevent=EEG.event(source_ev_id).urevent;
%                     EEG.event(n_ev+add_ev_ct).epoch=a;
%                     
%                     if ~EEG.epoch(a).eventlatency{b}
%                         %only update bin counter if the time locking event
%                         %falls in the current bin
%                         bin_ct(bin)=bin_ct(bin)+1;
%                     end
                end
            end
        else
            %if only one event per epoch, a string is used instead
            %of a cell array
            n_ep_type=1;
        end
        
        add_ct=0;
        for b=1:n_ep_type,
            for bin=1:n_bin,
                if n_ep_type>1,
                    if ismember(EEG.epoch(a).eventtype{b},bintype{bin})
                        %Add bin info to EEG.epoch
                        add_ct=add_ct+1;
                        add_ev_ct=add_ev_ct+1;
                        EEG.epoch(a).eventtype{add_ct+n_ep_type}=['bin' num2str(bin)];
                        %Duplicate other information about the event
                        if isfield(EEG.epoch(a),'eventduration')
                            %Some versions of EEGLAB apparently don't have the
                            %eventduration field or don't automatically
                            %add it
                            EEG.epoch(a).eventduration{add_ct+n_ep_type}=EEG.epoch(a).eventduration{b};
                        end
                        EEG.epoch(a).eventlatency{add_ct+n_ep_type}=EEG.epoch(a).eventlatency{b};
                        EEG.epoch(a).event(add_ct+n_ep_type)=n_ev+add_ev_ct;
                        EEG.epoch(a).eventurevent{add_ct+n_ep_type}=EEG.epoch(a).eventurevent{b};
                        
                        source_ev_id=EEG.epoch(a).event(b);
                        %Note, it looks like EEGLAB adds an event for a urevent everytime it occurs in an epoch
                        
                        %Add bin info to EEG.event
                        EEG.event(n_ev+add_ev_ct).type=['bin' num2str(bin)];
                        EEG.event(n_ev+add_ev_ct).latency=EEG.event(source_ev_id).latency;
                        if isfield(EEG.event(source_ev_id),'duration'),
                            EEG.event(n_ev+add_ev_ct).duration=EEG.event(source_ev_id).duration;
                        end
                        EEG.event(n_ev+add_ev_ct).urevent=EEG.event(source_ev_id).urevent;
                        EEG.event(n_ev+add_ev_ct).epoch=a;
                        
                        
                        if ~EEG.epoch(a).eventlatency{b}
                            %only update bin counter if the time locking event
                            %falls in the current bin
                            bin_ct(bin)=bin_ct(bin)+1;
                        end
                    end
                else
                    if ismember(EEG.epoch(a).eventtype,bintype{bin})
                        %convert non cell array fields of EEG.epoch(1) into
                        %cell arrays
                        
                        % Note: EEG.epoch(a).event stays a vector even when
                        % multiple events are present in an epoch
                        
                        if isfield(EEG.epoch(a),'eventduration')
                            %Older versions of EEGLAB apparently don't have the
                            %eventduration field
                            evdurflag=1;
                            if ~iscell(EEG.epoch(a).eventduration)
                                temp=EEG.epoch(a).eventduration; %EEG.eventduration
                                EEG.epoch(a).eventduration=cell(1,1);
                                EEG.epoch(a).eventduration{1}=temp;
                            end
                        else
                            evdurflag=0;
                        end
                        
                        if ~iscell(EEG.epoch(a).eventlatency)
                            temp=EEG.epoch(a).eventlatency; %EEG.eventlatency
                            EEG.epoch(a).eventlatency=cell(1,1);
                            EEG.epoch(a).eventlatency{1}=temp;
                        end
                        
                        if ~iscell(EEG.epoch(a).eventtype)
                            temp=EEG.epoch(a).eventtype; %EEG.eventtype
                            EEG.epoch(a).eventtype=cell(1,1);
                            EEG.epoch(a).eventtype{1}=temp;
                        end
                        
                        if ~iscell(EEG.epoch(a).eventurevent)
                            temp=EEG.epoch(a).eventurevent; %EEG.eventurevent
                            EEG.epoch(a).eventurevent=cell(1,1);
                            EEG.epoch(a).eventurevent{1}=temp;
                        end
                        
                        %Add bin info to EEG.epoch
                        add_ct=add_ct+1;
                        add_ev_ct=add_ev_ct+1;
                        EEG.epoch(a).eventtype{add_ct+n_ep_type}=['bin' num2str(bin)];
                        %Duplicate other information about event
                        if evdurflag,
                            EEG.epoch(a).eventduration{add_ct+n_ep_type}=EEG.epoch(a).eventduration{b};
                        end
                        EEG.epoch(a).eventlatency{add_ct+n_ep_type}=EEG.epoch(a).eventlatency{b};
                        EEG.epoch(a).event(add_ct+n_ep_type)=n_ev+add_ev_ct;
                        EEG.epoch(a).eventurevent{add_ct+n_ep_type}=EEG.epoch(a).eventurevent{b};
                        
                        source_ev_id=EEG.epoch(a).event(b);
                        %Note, it looks like EEGLAB adds an event for a urevent everytime it occurs in an epoch
                        
                        %Add bin info to EEG.event
                        EEG.event(n_ev+add_ev_ct).type=['bin' num2str(bin)];
                        EEG.event(n_ev+add_ev_ct).latency=EEG.event(source_ev_id).latency;
                        if isfield(EEG.event(source_ev_id),'duration')
                            EEG.event(n_ev+add_ev_ct).duration=EEG.event(source_ev_id).duration;
                        end
                        EEG.event(n_ev+add_ev_ct).urevent=EEG.event(source_ev_id).urevent;
                        EEG.event(n_ev+add_ev_ct).epoch=a;
                        
                        if ~EEG.epoch(a).eventlatency{b}
                            %only update bin counter if the time locking event
                            %falls in the current bin
                            bin_ct(bin)=bin_ct(bin)+1;
                        end
                    end
                end
            end
        end
    end
end

%% Output summary to command line and save EEG variable
if VERBLEVEL>1,
    fprintf('\nNumber of events per bin.\n');
    for a=1:n_bin,
        fprintf('Bin %d (%s): %d\n',a,EEG.bindesc{a},bin_ct(a));
    end
end

if ~isempty(save_fname),
    if VERBLEVEL>0,
        fprintf('Saving EEG variable as the file: %s\n',save_fname);
    end
    [justpath, justname]=pathNname(save_fname);
    EEG=pop_saveset(EEG,'filename',justname,'filepath',justpath);
end
