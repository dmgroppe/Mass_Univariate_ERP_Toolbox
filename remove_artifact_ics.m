% remove_artifact_ics() - Removes the independent components (ICs) of an
%                         EEGLAB "EEG" variable that have been labeled as an
%                         artifact in EEG.ic_labels or EEG.reject.gcompreject.
%                         The artifacts are removed from EEG.data and
%                         EEG.data is optionally baselined.  Artifact
%                         status of labels in EEG.ic_labels is based on
%                         is_art.m.  Note that EEG.ic_labels is a field
%                         unique to Kutaslab data.
%
%
% Usage:
%  >> art_ics=remove_artifact_ics(bsln_wind,verblevel);
%
% Required Global Variable:
%  EEG - EEGLAB EEG variable.  EEG.data field will be modified by this
%        function.
%
% Optional Inputs:
%  bsln_wind - [start_time stop_time or NaN] Two element vector specifying 
%              baseline window in milliseconds or NaN.  If NaN, data are 
%              not baselined. Otherwise, the mean EEG at each channel and
%              epoch in this window will be removed from each epoch.
%  verblevel - an integer specifiying the amount of information you want
%              functions to provide about what they are doing during runtime.
%                Options are:
%                  0 - quiet, only show errors, warnings, and EEGLAB reports
%                  1 - stuff anyone should probably know
%                  2 - stuff you should know the first time you start working
%                      with a data set {default value if not already globally
%                      specified}
%                  3 - stuff that might help you debug (show all
%                      reports)
%
% Outputs:
%  art_ics - [integer vector] Indices of ICs labeled as artifacts
%
%
% Author:
% David Groppe
% Kutaslab, 9/2010

% Changes
% 8/23/2012 - NaN now a possible value for bsln_wind to avoid any
% baselining

function art_ics=remove_artifact_ics(bsln_wind,verblevel)

global EEG;
global VERBLEVEL;

if isempty(EEG.icawinv),
    %ICA hasn't been applied to these data
    art_ics=[];
    return
else
    
    if nargin<1
        bsln_wind=[];
    end
    
    if nargin<2,
        if ~isnan(bsln_wind)
            fprintf('NOT baselining data.\n')
        elseif (length(bsln_wind)~=2) || ~isnumeric(bsln_wind),
            if ~isempty(bsln_wind)
                error('Baseline window needs to be a two element vector composed of start and stop times (in ms).');
            end
        end
        
        if isempty(VERBLEVEL)
            VERBLEVEL=2;
        end
    else
        VERBLEVEL=verblevel;
    end
    
    if nargin>2,
        error('remove_artifact_ics.m accepts only two arguments.');
    end
    
    
    [n_chans, n_pts, n_epochs]=size(EEG.data);
    
    if ~isempty(bsln_wind),
        bsln_pts(1)=find_tpt(bsln_wind(1),EEG.times);
        bsln_pts(2)=find_tpt(bsln_wind(2),EEG.times);
    end
    
    fldnames=fieldnames(EEG);
    iclabels_used=0;
    for fn=1:length(fldnames),
        if strcmpi(fldnames{fn},'iclabels')
            iclabels_used=1;
        end
    end
    
    n_ics=size(EEG.icawinv,2);
    if iclabels_used,
        n=length(EEG.iclabels);
    else
        n=0;
    end
    art_ics=[];
    if VERBLEVEL>=2,
        fprintf('Removing artifact ICs from EEG set file: %s\n',EEG.setname);
        if n
            %n is only non-zero if there are IC labels
            disp('Artifact ICs:');
        end
    end
    for dg=1:n,
        if is_art(EEG.iclabels{dg}),
            art_ics=[art_ics dg];
            if VERBLEVEL>=2,
                fprintf('IC %d: %s\n',dg,EEG.iclabels{dg});
            end
        end
    end
    
    %Check EEG.iclabels artifact ICs against EEG.reject
    eeglab_rej=find(EEG.reject.gcompreject==1);
    if iclabels_used,
        dif1=setdiff(eeglab_rej,art_ics);
        dif2=setdiff(art_ics,eeglab_rej);
        if ~isempty(dif1) || ~isempty(dif2)
            if isempty(eeglab_rej),
                rej_ic_str='None';
            else
                rej_ic_str=int2str(eeglab_rej);
            end
            msg=['The ICs marked for rejection by EEGLAB in EEG.reject.gcompreject differ from those labeled as artifacts in EEG.iclabels.' ...);
                10 'Only the artifact labels in EEG.iclabels will be used.' 10 ...
                'The ICs labeled by artifacts by EEG.reject.gcompreject are: ' rej_ic_str];
            watchit(msg);
        end
    else
        art_ics=eeglab_rej;
    end
    
    if VERBLEVEL>=2,
        fprintf('%d total artifact ICs will be removed.\n',length(art_ics));
    end
    
    nonart_ics=setdiff(1:n_ics,art_ics);
    if ~isempty(nonart_ics),
        %zero IC artifact ICS
        unmix=EEG.icaweights*EEG.icasphere;
        fltr=EEG.icawinv(:,nonart_ics)*unmix(nonart_ics,:);
        EEG.data=fltr*reshape(EEG.data,n_chans,n_pts*n_epochs);
        
        if isempty(bsln_wind) || sum(isnan(bsln_wind)),
            EEG.data=reshape(EEG.data,n_chans,n_pts,n_epochs);
        else
            %baseline data
            if VERBLEVEL>=2,
                fprintf('Baselining data by removing mean EEG between %d and %d ms (time points %d and %d).\n', ...
                    bsln_wind(1),bsln_wind(2),bsln_pts(1),bsln_pts(2));
            end
            
            EEG.data=rmbase(EEG.data,n_pts,bsln_pts(1):bsln_pts(2));
            EEG.data=reshape(EEG.data,n_chans,n_pts,n_epochs);
        end
    else
        error('All ICs have been labeled as artifacts.');
    end
    
end