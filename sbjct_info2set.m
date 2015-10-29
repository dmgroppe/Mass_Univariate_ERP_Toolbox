function sbjct_info2set(in_setfile,sbjct_infofile,sbjct_id,out_setfile,forcewrite,verblevel)
%
% sbjct_info2set() - Imports information about a subject into the epochs
%                     of an EEGLAB .set file from a text file.
%                     The new information (e.g., age, sex) is added to the 
%                     EEG.event and EEG.epoch fields of the EEG struct
%                     variable.
%                     
% Usage:
%  >> sbjct_info2set(in_setfile,sbjct_infofile,sbjct_id,out_setfile,forcewrite,verblevel)
%
%
% Required Inputs:
%   in_setfile        = [string] the name of the source EEGLAB .set file
%                        (include the file's pathname unless it is in the
%                        current working directory).
%
%   sbjct_infofile     = [string] a space or tab delimited text file
%                        containing additional numeric information
%                        about the subject who generated the data in
%                        'in_setfile'.  The first row of the file is a
%                        header line with a one word description of each
%                        column in the file. Each row below the headerline    
%                        corresponds to a different experimental subject.
%                        The first column of the file specifies the subject
%                        ID (it can be a number or a string). Additional columns
%                        specificy information about the subject (e.g.,
%                        age, sex).  This information can be numeric or a
%                        string.
%                        
%   sbjct_id           = a number or string that identifies the subject
%                        who produced the data in 'in_setfile.'  The contents
%                        of 'sbjct_id' should match the contents of exactly 
%                        one cell in the first column of the file specified by
%                        'sbjct_infofile' (matching is NOT case sensitive).
%
% Optional Inputs:
%   out_setfile     = [string] the name of the EEGLAB .set file that will
%                      be written to (include the file's pathname unless
%                      it is in the current working directory) {default:
%                      same as in_setfile}
%
%   forcewrite      = ['on' | 'off'] If 'on', the function will
%                      automatically overwrite any pre-existing event
%                      information in the .set file that contradicts the
%                      event information in sbjct_infofile.  If 'off,' the
%                      user will be asked whether or not to overwrite any
%                      such pre-existing information. {default: 'off'}
%
%   verblevel       = an integer specifiying the amount of information you
%                      want functions to provide about what they are doing
%                      during runtime.
%                        Options are:
%                         0 - quiet, only show errors, warnings, and EEGLAB
%                             reports
%                         1 - stuff anyone should probably know
%                         2 - stuff you should know the first time you start
%                             working with a data set {default value if
%                             verblevel not already specified}
%                         3 - stuff that might help you debug (show all
%                             reports)
%
% Outputs:
%   Function outputs nothing in Matlab.  It simply writes/overwites a .set
%   file to disk.
%
%
% Example:
% >> in_fname='/homes/dgroppe/SANDBOX/PRSENT/prsent41.set';
% >> sbjct_infofile='/homes/dgroppe/SANDBOX/PRSENT/prsent_sbjctinfo.txt';
% >> evcode_info2set(in_fname,sbjct_infofile,'temp.set');
%
%
% Additional Notes:
%
% Don't use parenthesis in column header names.  Matlab interprets
% the name as a function call.
%
% Use NaN to fill cells for codes that don't have a value for a
% particular column
%
% Author:
% David Groppe
% Kutaslab, 10/2009
%

global VERBLEVEL

%Check Inputs
if nargin<3,
    error('Function sbjct_info2set.m requires at least three arguments.\n');
end

if ~ischar(in_setfile),
    error('Argument in_setfile needs to be a string that specifies an EEGLAB .set file.');
end

if ~ischar(sbjct_infofile),
    error('Argument sbjct_infofile needs to be a string that specifies a text file.');
end

if isnumeric(sbjct_id),
    sbjct_id=num2str(sbjct_id);
end

if nargin<4,
    out_setfile=in_setfile; %default: output file same as input file
    VerbReport('New set file will overwrite old set file (default behavior).', ...
        2, VERBLEVEL);
else
    if ~ischar(out_setfile),
        error('Argument out_setfile needs to be a string that specifies an EEGLAB .set file.');
    end
end

if nargin<5,
    forcewrite='off';
elseif (~strcmpi(forcewrite,'on') && ~strcmpi(forcewrite,'off'))
    error('Argument "forcewrite" needs to be set to ''on'' or ''off''.');
end

if nargin<6
    if isempty(VERBLEVEL),
        VERBLEVEL=2; %default
    end
else
    VERBLEVEL=verblevel;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD .set FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG=pop_loadset(in_setfile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD INFORMATION ABOUT EACH CLASS OF EXPERIMENTAL EVENT (e.g. TARGETS vs. STANDARDS) FROM A TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VerbReport(sprintf('Getting additional information about subject from %s', ...
    sbjct_infofile), 1, VERBLEVEL);
[sub_fid, message]=fopen(sbjct_infofile,'r');
if (sub_fid==-1),
    fprintf('*************** ERROR ******************\n');
    fprintf('Cannot open file %s.\n',sbjct_infofile);
    fprintf('According to fopen: %s.\n',message);
    error('Aborting import of additional subject information.');
else
    %read column headers
    txtline = fgetl(sub_fid);
    if (txtline==-1),
        fprintf('*************** ERROR ******************\n');
        fprintf('File %s is empty.\n',sub_fid);
        error('Aborting import of additional subject information.');
    else
        
        %Parse column header
        clear sub_col_hdrs;
        col_ct=1;
        [sub_col_hdrs{col_ct}, rmndr]=strtok(txtline);
        fprintf('Subject ID column is: %s\n',sub_col_hdrs{col_ct});
        while ~isempty(rmndr) && charleft(rmndr),
            col_ct=col_ct+1;
            [sub_col_hdrs{col_ct}, rmndr]=strtok(rmndr);
            fprintf('Column %d is: %s\n',col_ct,sub_col_hdrs{col_ct});
        end
        n_col=col_ct;
        
        %Read rows (information per subject)
        row_ct=1;
        found_flag=0;
        while ~feof(sub_fid)
            txtline = fgetl(sub_fid);
            [neo_val, txtline]=strtok(txtline);
            if strcmpi(neo_val,sbjct_id),
                if found_flag,
                    error('File %s had multiple rows for subject %s.  It should only have one.',sbjct_infofile,sbjct_id);
                else
                    found_flag=1;
                    col_ct=1;
                    file_sub_info=cell(1,n_col-1);
                    while ~isempty(txtline) && charleft(txtline),
                        [neo_val, txtline]=strtok(txtline);
                        file_sub_info{col_ct}=neo_val; %novel info is a string
                        %Note, subjects that do not have a value for that
                        %column should be represented as NaN
                        col_ct=col_ct+1;
                    end
                end
            end
            row_ct=row_ct+1;
        end
    end
    fclose(sub_fid);
end
if ~found_flag,
   error('Subject %s not found in file %s.  Import of subject information aborted.', ...
       sbjct_id,sbjct_infofile);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD NEW EPOCH INFORMATION TO EEG STRUCT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fldnames=fieldnames(EEG.epoch);
%for each new column (assume first column is Subject ID)
for a=2:length(sub_col_hdrs),
    doit=1;
    %check to make sure column does not already exist
    for b=1:length(fldnames),
        if strcmpi(['event' sub_col_hdrs{a}],fldnames{b}),
            if ~strcmpi(forcewrite,'on')
                %ask user if s/he wants to overwrite
                fprintf('Field %s already exists in %s.\n',sub_col_hdrs{a},EEG.filename);
                erayz=input('Do you want to overwrite existing field? ("y" or "n"): ','s');
                while ~strcmpi(erayz,'y') && ~strcmpi(erayz,'n')
                    erayz=input('Please enter the single letter "y" or "n": ','s');
                end
            else
                erayz='y';
            end
            if strcmpi(erayz,'n')
                fprintf('NOT overwriting EEG.%s.\n',fldnames{b});
                doit=0;
                break;
            else
                fprintf('OVERWRITING EEG.%s.\n',fldnames{b});
                EEG.epoch=rmfield(EEG.epoch,fldnames{b});
                EEG.event=rmfield(EEG.event,fldnames{b}(6:end)); %starting at 6 ignores 'event'
                break;
            end
        end
    end
    
    %write new event info to struct
    if doit,
        if isempty(str2num(file_sub_info{a-1})),
            str=1; %string
        else
            str=0; %numeric
        end
        
        %add to EEG.epoch
        for c=1:length(EEG.epoch),
            if str,
                cmd=['EEG.epoch(' int2str(c) ').event'  sub_col_hdrs{a} '='''  file_sub_info{a-1} ''';']; 
            else
                cmd=['EEG.epoch(' int2str(c) ').event'  sub_col_hdrs{a} '='  file_sub_info{a-1} ';']; 
            end
            eval(cmd);
        end
        
        %add to EEG.event
        for c=1:length(EEG.event),
            if str,
                cmd=['EEG.event(' int2str(c) ').'  sub_col_hdrs{a} '='''  file_sub_info{a-1} ''';'];%string
            else
                cmd=['EEG.event(' int2str(c) ').'  sub_col_hdrs{a} '='  file_sub_info{a-1} ';'];%numeric
            end
            eval(cmd);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE .set FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fpath, fname]=pathNname(out_setfile);
EEG=pop_saveset(EEG,'filepath',fpath,'filename',fname);

end


function yesno=charleft(str)


yesno=0;
for x=1:length(str),
    if (str(x)~=32) && (str(x)~=13)
        %32=space, 13=carriage return
        yesno=1;
        break;
    end
end


end