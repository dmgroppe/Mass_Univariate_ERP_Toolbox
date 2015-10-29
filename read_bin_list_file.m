% read_bin_list_file() - Extracts "bin" information from a text file called
%         a "bin list file."  The bin list file (blf) indicates which types 
%         of events belong to which bins.  Each line of the blf needs to 
%         start with a number followed by a closed parenthesis, which 
%         indicates the bin number.  Next comes a string, number, or a 
%         series of numbers indicting the "type" of events that fall in 
%         that bin (note the "type" of an event is stored in the field 
%         EEG.epoch(#).eventtype and EEG.event(#).type of an EEGLAB EEG 
%         variable).  Finally, comes an equal sign and text that describes
%         the kinds of trials that fall into that bin.
%
% For example, the following line would classify epochs of type 22 as
% members of Bin 1, which is described as containing "Targets in Auditory Task".
% 1) 22=Targets in Auditory Task
%
% The following would classify epochs of type 22 and 26 as members of Bin 1, 
% which is described as containing "Targets in Auditory Task".
% 1) 22 26=Targets in Auditory Task
%
% The following would classify epochs of type 22, 23, 24 and 25 as members 
% of Bin 1, which is described as containing "Targets in Auditory Task" 
% (note that using MATLAB notation to indicate vectors is acceptable).
% 1) 22:25=Targets in Auditory Task
%
% The following would classify epochs of type 'aud_targ' as members of Bin 1, 
% which is described as containing "Targets in Auditory Task".
% 1) aud_targ=Targets in Auditory Task
%
% The following would classify epochs of type 'aud targ' as members of Bin 1, 
% which is described as containing "Targets in Auditory Task". Note, the
% use of single quotes to include the space in the type name.  Without the
% single quotes, two types would be read: "aud" and "targ".
% 1) 'aud targ'=Targets in Auditory Task
%
%
% Usage:
%  >> [bindesc, bintype]=read_bin_list_file(blf_fname, verblevel);
%
%
% Required Input:
%   blf_fname  - [string] A set of integers specifying the subject ID 
%                numbers to include in the grand average.  Only necessary 
%                if a filename template is given as the input to 
%                gui_infiles_or_tmplt.
%
% Optional Input:
%   verblevel        - An integer specifiying the amount of information you want
%                      this function to provide about what it is doing during runtime.
%                       Options are:
%                        0 - quiet, only show errors, warnings, and EEGLAB reports
%                        1 - stuff anyone should probably know
%                        2 - stuff you should know the first time you start working
%                            with a data set {default value if not globally specified}
%                        3 - stuff that might help you debug (show all
%                            reports)
%
% Output:
%   bindesc  - A cell array of strings describing the contents of each bin.
%   bintype  - A cell array of cell arrays of strings listing the types of
%              events that fall in each bin
%                
% Global Variable:
%   VERBLEVEL - Mass Univariate ERP Toolbox level of verbosity (i.e., tells 
%               functions how much to report about what they're doing during
%               runtime) set by the optional function argument 'verblevel'
%
%
% NOTES:
% -You can use single quotes to define type names that have spaces in them
% (e.g., 'S  1').  However you may only do this once per line of your blf
% file.  In other words you can have exactly zero or two single quotes per
% line.
%
% Author:
% David Groppe
% Kutaslab, 9/2010

%%%%%%%%%%%%%%%% Future Work  %%%%%%%%%%%%%%%%%
%

function [bindesc, bintype]=read_bin_list_file(blf_fname, verblevel)

if nargin<2,
   verblevel=2; 
end

[fid, msg]=fopen(blf_fname);
if fid==-1,
    error('Couldn''t open blf file %s because: %s',blf_fname,msg);
end
line_ct=0;
n_bin=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline),
        break; %leave while loop
    end
    line_ct=line_ct+1;
    [bin_str tline]=strtok(tline);
    bin_id(line_ct)=str2double(bin_str(1:length(bin_str)-1)); %str2double returns NaN if the string it's converting contains alpha numeric characters
    if isnan(bin_id(line_ct)),
        error('Line #%d of file %s should start with a bin number followed by a close parenthesis\n',line_ct, ...
            blf_fname,line_ct);
    end
    eql_id=0;
    for a=1:length(tline),
        if tline(a)=='=',
            eql_id=a;
            break;
        end
    end
    if eql_id,
        line_type=rm_leading_whitespace(tline(1:eql_id-1));
        line_type_as_num=str2num(line_type); %str2num returns the empty set 
        %if the string cannot be interpreted as numbers at the MATLAB command line (e.g., ':')
        if ~isempty(line_type_as_num),
           line_type=num2str(line_type_as_num); %convert it back to string 
        end
        if line_type(1)==39,
           if sum(line_type==39)~=2,
              error('Problem with Line %d. You can only use a pair of single quotes to indicate a type.',line_ct);  
           end
           quote_locs=find(line_type==39);
           type_array{1}=line_type(2:quote_locs(2)-1);
        else
            type_array=strread(line_type,'%s');
        end
    else
        error('Line %d needs to have an equals sign\n',line_ct);
    end
    
    if (bin_id(line_ct)>n_bin) || isempty(bindesc{bin_id(line_ct)})
        %new bin
        bindesc{bin_id(line_ct)}=rm_leading_whitespace(tline(eql_id+1:end));
        n_bin=n_bin+1;
        n_type=length(type_array);
        for x=1:n_type,
           bintype{bin_id(line_ct)}{x}=type_array{x}; 
        end
    else
       %The bin for this line has already been entered, make sure it's
       %consistent
       neo_bindesc=rm_leading_whitespace(tline(eql_id+1:end));  
       if ~strcmpi(neo_bindesc,bindesc{bin_id(line_ct)})
           error('Bin descriptor on Line %d (%s) does not equal previous descriptor for that bin (i.e., Bin %d: %s)', ...
               line_ct,neo_bindesc,bin_id(line_ct),bindesc{bin_id(line_ct)});
       end
       n_old_type=length(bintype{bin_id(line_ct)});
       n_new_type=length(type_array);
       for x=1:n_new_type,
           bintype{bin_id(line_ct)}{x+n_old_type}=type_array{x};
       end
    end
end
fclose(fid);

%% Error checking

%Make sure Bin ID #'s are contiguous and start with 1
uni_bin=unique(bin_id);
if max(abs(uni_bin-[1:n_bin])),
   error('Bin numbers need to be contiguous and start a "1".  Your bin numbers are: %s\n',num2str(uni_bin));
end

%Get rid of redundant types (just in case)
for a=1:n_bin,
    bintype{a}=unique(bintype{a});
end
     
if verblevel>1,
    %Print which event types belong to each bin for user to double check
    for a=1:n_bin,
        fprintf('\nBin %d: %s\n',a,bindesc{a});
        fprintf('Event types belonging to this bin:')
        for b=1:length(bintype{a})
            fprintf(' ''%s''',bintype{a}{b});
        end
        fprintf('\n')
    end
end


function str=rm_leading_whitespace(str)

n=length(str);
for a=1:n,
   if str(a)~=' ',
      str=str(a:end);
      break;
   end
end