function reject=is_art(label)
% isart() - Function for evaluating whether or not the text label of an
%           independent component (IC) corresponds to an EEG artifact or not.
%              
% Usage:
%  >> reject=isart(label)
%
% Required Inputs:
%  label - a string label for an IC.
%
% Output:
%  reject - 1 if "label" is a type of IC artifact, 0 otherwise. The
%            possible values of label recognized as IC artifacts are:
%              'blink'    -Eye blink IC
%              'HE'       -Horizontal eye movement IC
%              'EOG'      -Generic eye movment IC
%              'EMG'      -Muscle noise IC
%              'heart'    -Heartbeat artifact IC
%              'drift'    -Drift artifact IC
%              'bad_chan' -Bad channel (i.e., sensor problems) IC
%              '60HZ'     -60 Hz line noise IC
%              'art'      -Generic artifact IC
%              'blink (objective)' -Eye blink IC (according to objective test)
%              'heart (objective)' -Heart blink IC (according to objective
%                                   test)
%              'HE (objective)'    -Horizontal eye movement IC (according to objective
%                                   test)
%              'EMG (objective)'   -Muscle noise IC (according to objective
%                                   test)
%              'bad_chan (objective)' -Bad channel (i.e., sensor problems)
%                                      IC (according to objective test)
%              'xtrm_activity (objective)' -Extreme activations (most likely 
%                                           due to movement of some sort)
%              'EOG (objective)' -Eye movement IC (according to objective
%                                 test)
%
% Notes:
% -The evaluation of "label" is NOT case sensitive.  For example,
% 'blink','Blink','BLINK', and 'bLiNK' would all be recognized as an
% artifact.
%
% Author: 
% David Groppe
% Kutaslab, 11/2009 
%


art_labels{1}='blink';
art_labels{2}='HE';
art_labels{3}='EOG';
art_labels{4}='NZ'; %left in for backwards compatibility
art_labels{5}='heart';
art_labels{6}='drift';
art_labels{7}='bad_chan';
art_labels{8}='60HZ';
art_labels{9}='art';
art_labels{10}='EMG';
art_labels{11}='blink (objective)';
art_labels{12}='heart (objective)';
art_labels{13}='HE (objective)';
art_labels{14}='EMG (objective)';
art_labels{15}='bad_chan (objective)';
art_labels{16}='xtrm_activity (objective)';
art_labels{17}='EOG (objective)';


reject=0;
if ~isempty(label)
    for a=1:length(art_labels),
        if strcmpi(label,art_labels{a})
            reject=1;
            return;
        end
    end
end
