% list_event_types() - Tallies the set of event "types" in an EEGLAB EEG
%                      variable and the number of times each event occurs.
%
% Usage:
%  >>  [types type_count]=list_event_types(EEG_or_fname,timelock_only)
%
%
% Required Input:
%   EEG_or_fname - EEGLAB EEG struct variable containing epoched or 
%                  continuous EEG data or the filename of a set file
%                  that contains such an EEG struct variable.
%
%
% Optional Input:
%   timelock_only - [1 or 0] If 1, only the event types of events to which
%                   an epoch is time locked will be tallied (i.e., events 
%                   that occur at 0 ms in the epoch). If 0, all events 
%                   types are tallied. This optional input is only relevant
%                   to epoched data sets. {default: 1}
%
% Outputs:
%   types       - Cell array of the set of event types found  
%
%   type_count  - Number of instances of each type found. Note that
%                 type_count(a) is the count for types{a}.
%
% Author:
% David Groppe
% Kutaslab, 10/2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 12/16/2010-Can now handle continuous data
%
% 7/10/2011-Can now deal with event types that aren't strings

function [types type_count]=list_event_types(EEG_or_fname,timelock_only)

if nargin<2,
   timelock_only=1; 
end

if ischar(EEG_or_fname),
    EEG=pop_loadset(EEG_or_fname);
else
    EEG=EEG_or_fname;
    clear EEG_or_fname;
end

%% Modicum of Error Checking
if isempty(EEG.epoch),
    %Continuous data
    fprintf('set file contains continuous EEG data.\n');
    n_ev=length(EEG.event);
    n_types=0;
    types=[];
    type_count=[];
    for a=1:n_ev,
        ev_type=EEG.event(a).type;
        if ~isstr(ev_type)
            ev_type=num2str(ev_type);
        end
        if isempty(types),
            tf=0;
        else
            [tf, loc]=ismember(ev_type,types);
        end
        if tf
            %type has already been encountered, increment counter
            type_count(loc)=type_count(loc)+1;
        else
            %novel type
            n_types=n_types+1;
            types{n_types}=ev_type;
            type_count(n_types)=1;
        end
    end
else
    %Epoched data
    fprintf('set file contains epoched EEG data.\n');
    n_ep=length(EEG.epoch);
    
    % Check for empty epochs. Remove them if found.
    non_empty_ids=[];
    for a=1:n_ep,
        if ~isempty(EEG.epoch(a).event)
            non_empty_ids=[non_empty_ids a];
        end
    end
    n_ep_orig=n_ep;
    EEG.epoch=EEG.epoch(non_empty_ids);
    n_ep=length(EEG.epoch);
    if n_ep_orig>n_ep,
        warning(sprintf('%d empty epoch(s) removed',n_ep_orig-n_ep));
    end
    
    n_types=0;
    types=[];
    type_count=[];
    if timelock_only,
        for a=1:n_ep,
            if iscell(EEG.epoch(a).eventtype),
                %if only one event per epoch, a string is used instead
                %of a cell array
                n_type_this_ep=length(EEG.epoch(a).eventtype);
            else
                n_type_this_ep=1;
            end
            for b=1:n_type_this_ep,
                if iscell(EEG.epoch(a).eventlatency),
                    %if only one event per epoch, a scalar is used instead
                    %of a cell array
                    ev_latency=EEG.epoch(a).eventlatency{b};
                    ev_type=EEG.epoch(a).eventtype{b};
                else
                    ev_latency=EEG.epoch(a).eventlatency(b);
                    ev_type=EEG.epoch(a).eventtype;
                end
                
                if ~ev_latency,
                    %only log event types for events that the epoch is
                    %timelocked to
                    if isempty(types)
                        tf=0;
                    else
                        [tf loc]=ismember(ev_type,types);
                    end
                    if tf
                        %type has already been encountered, increment counter
                        type_count(loc)=type_count(loc)+1;
                    else
                        %novel type
                        n_types=n_types+1;
                        if iscell(ev_type),
                            types{n_types}=ev_type{1};
                        else
                            types{n_types}=ev_type;
                        end
                        type_count(n_types)=1;
                    end
                end
            end
        end
    else
        for a=1:n_ep,
            if iscell(EEG.epoch(a).eventtype),
                %if only one event per epoch, a string is used instead
                %of a cell array
                n_type_this_ep=length(EEG.epoch(a).eventtype);
            else
                n_type_this_ep=1;
            end
            for b=1:n_type_this_ep,
                if iscell(EEG.epoch(a).eventlatency),
                    %if only one event per epoch, a string is used instead
                    %of a cell array
                    ev_type=EEG.epoch(a).eventtype{b};
                else
                    ev_type=EEG.epoch(a).eventtype;
                end
                [tf loc]=ismember(ev_type,types);
                if tf
                    %type has already been encountered, increment counter
                    type_count(loc)=type_count(loc)+1;
                else
                    %novel type
                    n_types=n_types+1;
                    types{n_types}=ev_type;
                    type_count(n_types)=1;
                end
            end
        end
    end
end

for a=1:n_types,
    fprintf('Event Type: ''%s'' (occurences=%d)\n',types{a},type_count(a));
end
