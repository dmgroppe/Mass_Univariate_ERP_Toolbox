% GNDs2GRP() - Creates a Mass Univariate ERP Toolbox GRP variable from a set
%              of Mass Univariate ERP Toolbox GND variables.  A GRP variable 
%              is useful for doing between-subjects analyses.  Note, all the 
%              GND files need to have samples at the same points in time, 
%              the same baseline window, and compatible sets of electrodes.
%             
% Usage:
%  >> GRP=GNDs2GRP(gui_infiles_or_tmplt,varargin);
%
% Required Inputs:
%   gui_infiles_or_tmplt - ['gui', a cell array of strings, or a string template]
%                          If 'gui', a GUI is created that allows you to select
%                          which average files to import (this is probably the
%                          easiest way to import files).  If a cell array of 
%                          of strings, each element of the cell array should 
%                          contain the filename of a GND file (e.g., 
%                          {'visodbl_young.GND','visodbl_old.GND'}.  Otherwise, this 
%                          input should be a filename template (i.e., a string with
%                          a # where an ID number should be--'visodbl#.GND').  
%                          If you provide a template, you must use the option
%                          'GND_ids' (see below) to specify the ID numbers of
%                          each GND file. Include the files' path unless the files 
%                          are in the current working directory.
%
% Optional Inputs:
%   group_desc       - [cell array] A string for each GND file that
%                      describes that group of participants (e.g., {'young',
%                      'old'}).  It might save you typing later to keep 
%                      these names simple.  If this option isn't used, a 
%                      pop-up window will be created for each GND file to 
%                      request a name for that group. {default: not used}
%
%   exp_desc         - [string] A string that described the experiment. If
%                      not specified, the experiment descriptor from the 
%                      first GND file is used. {default: not specified}
%
%   GND_ids          - [integer vector] A set of integers specifying the
%                      GND file ID numbers of each group of participants.  
%                      Only necessary if a filename template is given as 
%                      the input to gui_infiles_or_tmplt. {default: not
%                      used}
%
%   create_difs      - ['yes' or 'no'] If 'yes', difference waves will be
%                      created for each bin in the GND files between all 
%                      possible pairs of groups and added as new bins to 
%                      the GRP variable. The GND files must all have the 
%                      exact same bins (including bin descriptors) to use 
%                      this option. {default: not used}
%
%   use_bins         - [integer vector] A set of integers specifying which
%                      bins to take difference waves.  If not specified, all 
%                      bins will be imported. Note, if you import only a 
%                      subset of bins, the bin numberings will start at 1 
%                      and run to the number you've imported (i.e., the may
%                      differ from the bin numbers in the average file. 
%                      {default: import all bins}
%
%   exclude_chans    - A cell array of channel labels to exclude from 
%                      importing (e.g., {'A2','lle','rhe'}). You cannot
%                      use both this option and 'include_chans' (below).{default:
%                      not used, import all channels}
%
%   include_chans    - A cell array of the labels of the channels you wish
%                      to import (e.g., {'A2','lle','rhe'}).  All other 
%                      channels will be ignored. You cannot use both this  
%                      option and 'exclude_chans' (above). {default: not  
%                      used, import all channels}
%
%   out_fname        - [string] Filename to save GND variable to.  If empty
%                      (i.e., not specified), a GUI will be created to prompt
%                      you for a filename. If the string 'no save', the GND 
%                      variable will not be saved to disk. 
%
% Output:
%   GRP  - Struct variable containing grand averages derived from pairs of 
%          GND variables, the file names of locations of the source GND files, 
%          t-test results and more.
%
% Example:
%   See documentation on between-participant analyses:
%      http://openwetware.org/wiki/Mass_Univariate_ERP_Toolbox:_between-subject_t-tests
%   or
%      http://kutaslab.pbworks.com/Between_Subject_Perm_Test
%
% Notes:
% -Individual participant ERPs are not stored with the GRP variable.
% Rather, whenever they are needed by MATLAB, the source GND variables are
% temporarily reloaded.  Thus, it is critical for operations like the
% creation of new bins, that the source GND variable files keep the same
% file name and location as when the GRP was first created.
%
% Author:
% David Groppe
% Kutaslab, 3/2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 5/7/2010-GRP.perm_tests turned into GRP.t_tests to be compatible with FDR
% code
%
% 10/25/2010-Modified to be able to deal with GND variables with no cal
% pulse information.  If a single GND variable does not have cal pulse
% information, the GRP variable will not have cal pulse information
% (GRP.cals will be empty)
%
% 4/4/2011: Made compatible with Windows.
%

%%%%%%%%%%%%%%%% FUTURE WORK %%%%%%%%%%%%%%%%%
% -Add Verblevel optional input?
% -Would it be useful to add the field GRP.sub_ct?
% -Make sure one can write Kutaslab avg files from a GRP variable
% -Is it really useful to have cal pulse information in a GRP variable?
% -Would it be useful to have the field GRP.indiv_bin_erps?

function GRP=GNDs2GRP(gui_infiles_or_tmplt,varargin)

p=inputParser;
p.addRequired('gui_infiles_or_tmplt',@(x) ischar(x) || iscell(x));
p.addParamValue('GND_ids',[],@isnumeric);
p.addParamValue('use_bins',[],@isnumeric);
p.addParamValue('exclude_chans',[],@(x) ischar(x) || iscell(x)); 
p.addParamValue('include_chans',[],@(x) ischar(x) || iscell(x));
p.addParamValue('out_fname',[],@ischar);
p.addParamValue('exp_desc',[],@ischar);
p.addParamValue('group_desc',[],@iscell);
p.addParamValue('create_difs','no',@(x) ischar(x) && (strcmpi(x,'no') || strcmpi(x,'yes')) );
p.parse(gui_infiles_or_tmplt,varargin{:});


%Figure out which channels to ignore if any
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

%Initialize GRP
GRP.exp_desc=p.Results.exp_desc;
GRP.filename=[];
GRP.filepath=[];
GRP.saved='no';
GRP.group_desc=[];
GRP.GND_fnames=[];
GRP.grands=[];   % <-grand average ERPs (chan x time point x  bin)
GRP.grands_stder=[]; % <-standard error of grand average ERPs (chan x time point x  bin)
GRP.grands_t=[]; %<-grand average ERPs as t-scores (chan x time point x  bin)
GRP.chanlocs=[];
GRP.bin_info=[]; %<-stores information about how difference waves were derived (which groups and bins)
GRP.condesc=[];
GRP.time_pts=[];
GRP.bsln_wind=[];
GRP.odelay=[];
GRP.srate=[];
GRP.indiv_fnames=[];
GRP.indiv_subnames=[];
GRP.indiv_traits=[];
GRP.indiv_bin_ct=[];
GRP.indiv_bin_raw_ct=[];
GRP.indiv_art_ics=[];
GRP.cals.indiv_cals=[]; 
GRP.cals.indiv_cal_ct=[];
GRP.cals.grand_cals=[];
GRP.cals.caldesc=[]; 
GRP.cals.condcode=0;
GRP.cals.condesc=[];
GRP.history=[];
GRP.t_tests=[];


if strcmpi(gui_infiles_or_tmplt,'GUI'),
    loading=1;
    infiles=[];
    while loading,
        fprintf('Select GND files to import.\n');
        fprintf('To select multiple contiguous files in the same directory, hold down the shift key when selecting.\n');
        fprintf('To select files in multiple directories or non-contiguous files in the same direcotry,\n');
        fprintf('select all the files in one directory and hit the "Open" button.\n');
        fprintf('You will then be given the opportunity to load more files from the current or another directory.\n');
        [neofname, inpath]=uigetfile({'*.GND','*.GND files'; ...
            '*.*','All Files (*.*)'},'Mass Uni GND Files to Import','MultiSelect','on');
        if ischar(neofname),
            clear infname;
            infname{1}=neofname; %make it a cell array for consistent syntax below
        else
            infname=neofname;
        end
        if ~inpath,
            if isempty(infiles),
                fprintf('File selection cancelled.  Aborting GNDs2GRP.\n');
                GRP=[];
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
    if isempty(p.Results.GND_ids)
        error('If gui_infiles_or_tmplt is a string filename template, you must specify the GND file numbers with the argument ''GND_ids''.');
    elseif ~ismember('#',p.Results.gui_infiles_or_tmplt),
        error('The filename template %s needs to contain a # to indicate where the GND file numbers go (e.g., yngvob#.GND).', ...
            p.Results.gui_infiles_or_tmplt);
    else
        lb_id=find(p.Results.gui_infiles_or_tmplt=='#');
        prefix=p.Results.gui_infiles_or_tmplt(1:(lb_id(1)-1)); %indexing lb_id by 1 in case there are multiple pound signs
        postfix=p.Results.gui_infiles_or_tmplt((lb_id(1)+1):end);
        n_subs=length(p.Results.GND_ids);
        infiles=cell(1,n_subs);
        for s=1:n_subs,
            no_pad=[prefix num2str(p.Results.GND_ids(s)) postfix];
            if p.Results.GND_ids(s)<10,
                padded=[prefix num2str(p.Results.GND_ids(s),'%.2d') postfix];
                if ismac || isunix
                    %Mac or Unix OS
                    [sP wP]=unix(['ls ' padded]); %if sp==0, then the file exists, need wP argument so that it isn't automatically displayed in command window
                    [sNP wNP]=unix(['ls ' no_pad]);
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

% Add group descriptors if specified
if ~isempty(p.Results.group_desc),
    if length(p.Results.group_desc)~=length(infiles),
        error('The number of group descriptors needs to equal the number of GND files.');
    else
        GRP.group_desc=p.Results.group_desc;
    end
end

if strcmpi(p.Results.create_difs,'yes'),
    create_difs=1;
else
    create_difs=0;
end

GRP.GND_fnames=infiles;
n_groups=length(infiles);
if n_groups<2,
   error('You need at least two GND files to make a GRP variable.'); 
end

for grp=1:n_groups,
    
    %load new GND variable
    GND=[];
    load(infiles{grp},'-MAT');
    if isempty(GND),
       error('File %s needs to contain a Mass Univariate ERP Toolbox GND variable that is called "GND"',infiles{grp}); 
    end
       
    if isempty(p.Results.group_desc),
       resp=inputdlg(sprintf('Enter the group descriptor for the participants in file %s (e.g., patients).', ...
           infiles{grp}),'Name Group of Participants');
       GRP.group_desc{grp}=resp{1};
    end
    
    %Set group general fields/check for GND field compatibility
    if grp==1,
        %set group general fields
        [n_chans, n_tpts, n_bins, n_subs1]=size(GND.indiv_erps);
        
        if create_difs,
            % select bins (in case only a subset were requested)
            if isempty(p.Results.use_bins),
                use_bins=1:n_bins;
            else
                use_bins=p.Results.use_bins;
            end
            n_use_bins=length(use_bins);
        end
            
        if isempty(GRP.exp_desc),
            GRP.exp_desc=GND.exp_desc;
        end        
        
        % select channels (in case only a subset were requested)
        if ~isempty(exclude_chans),
            use_chans=ones(1,n_chans); %preallocate mem
            for c=1:n_chans,
                for a=1:length(exclude_chans)
                    if strcmpi(GND.chanlocs(c).labels,exclude_chans{a}),
                        use_chans(c)=0;
                        break;
                    end
                end
            end
            use_chans=find(use_chans>0);
        elseif ~isempty(include_chans),
            use_chans=zeros(1,n_chans); %preallocate mem
            for c=1:n_chans,
                for a=1:length(include_chans)
                    if strcmpi(GND.chanlocs(c).labels,include_chans{a}),
                        use_chans(c)=1;
                        break;
                    end
                end
            end
            use_chans=find(use_chans>0);
        else
            use_chans=1:n_chans;
        end
        GRP.chanlocs=GND.chanlocs(use_chans);
        
        if create_difs,
            GRP.condesc=GND.condesc;
        end
        GRP.time_pts=GND.time_pts;
        GRP.bsln_wind=GND.bsln_wind;
        GRP.odelay=GND.odelay;
        GRP.srate=GND.srate;
        if isempty(GND.cals),
            GRP.cals=[];
        else
            GRP.cals.grand_cals=zeros(length(use_chans),n_tpts,n_groups);
            GRP.cals.caldesc=GND.cals.caldesc;
            GRP.cals.condesc=GND.cals.condesc;
        end
        
        if create_difs,
            %preallocate mem
            bin_sub_ct=zeros(n_use_bins,n_groups);
            bindesc=cell(1,n_use_bins);
            bin_condcodes=zeros(1,n_use_bins);
            bin_ct=0;
            for b=use_bins,
                bin_ct=bin_ct+1;
                bindesc{bin_ct}=GND.bin_info(b).bindesc;
                bin_condcodes(bin_ct)=GND.bin_info(b).condcode;
            end
        end
    else
        %check for compatibility with previously loaded GND variables
        [n_chansNEO, n_tptsNEO, n_binsNEO, n_subsNEO]=size(GND.indiv_erps);
                
        % select channels (in case only a subset were requested)
        if ~isempty(exclude_chans),
            use_chans=ones(1,n_chansNEO); %preallocate mem
            for c=1:n_chansNEO,
                for a=1:length(exclude_chans),
                    if strcmpi(GND.chanlocs(c).labels,exclude_chans{a}),
                        use_chans(c)=0;
                        break;
                    end
                end
            end
            use_chans=find(use_chans>0);
        elseif ~isempty(include_chans),
            use_chans=zeros(1,n_chansNEO); %preallocate mem
            for c=1:n_chansNEO,
                for a=1:length(include_chans)
                    if strcmpi(GND.chanlocs(c).labels,include_chans{a}),
                        use_chans(c)=1;
                        break;
                    end
                end
            end
            use_chans=find(use_chans>0);
        else
            use_chans=1:n_chansNEO;
        end
        
        %Make sure new GND variable in consistent with the old
        if ~isequal(GND.chanlocs(use_chans),GRP.chanlocs),
            error('The channel location information in the GND variable from file %s differs from that in previous files.', ...
                infiles{grp});
        end
        
        if create_difs,
            if ~isequal(GND.condesc,GRP.condesc),
                error('The condition code descriptor(s) in the GND variable from file %s differs from those in previous files.', ...
                    infiles{grp});
            end
        end
          
        if ~isequal(GND.time_pts,GRP.time_pts),
            error('The time points in the GND variable from file %s differ from those in previous files.', ...
                infiles{grp});
        end
        
        if ~isequal(GND.bsln_wind,GRP.bsln_wind),
            error('The baseline time window in the GND variable from file %s differs from that in previous files.', ...
                infiles{grp});
        end
        
        if ~isequal(GND.odelay,GRP.odelay),
            error('The odelay in the GND variable from file %s differs from that in previous files.', ...
                infiles{grp});
        end
        
        if ~isequal(GND.srate,GRP.srate),
            error('The sampling rate in the GND variable from file %s differs from that in previous files.', ...
                infiles{grp});
        end
        
        if ~isempty(GND.cals) && ~strcmpi(GND.cals.caldesc,GRP.cals.caldesc)
            watchit(sprintf('The cal pulse bin descriptor in the GND variable from file %s differs from that in previous files.', ...
                infiles{grp}));
        end
        
        if ~isempty(GND.cals) && ~strcmpi(GND.cals.condesc,GRP.cals.condesc)
            watchit(sprintf('The cal pulse condition code descriptor in the GND variable from file %s differs from that in previous files.', ...
                infiles{grp}));
        end
               
        if create_difs,
            bin_ct=0;
            for b=use_bins,
                bin_ct=bin_ct+1;
                if ~strcmpi(bindesc{bin_ct},GND.bin_info(b).bindesc),
                    watchit('Bin descriptor for Bin %d in file %s, differs from that in previous file(s).', ...
                        b,infiles{grp});
                    fprintf('Bin descriptor for Bin %d, file %s is %s.\n',b,infiles{grp},GND.bin_info(b).bindesc)
                    fprintf('Bin descriptor for Bin %d, previous file(s) is %s. Using bin descriptor from previous file(s).\n',b,bindesc{bin_ct});
                end
                if bin_condcodes(bin_ct)~=GND.bin_info(b).condcode,
                    watchit('Condition code for Bin %d in file %s is %d, but in previous file(s) it is %d.', ...
                        b,infiles{grp},GND.bin_info(b).condcode,bin_condcodes(bin_ct));
                    fprintf('Using condition code from previous file(s).\n')
                end
            end
        end
    end
    if ~isempty(GND.cals) && (GND.cals.condcode~=0)
        watchit(sprintf('The cal pulse condition code in the GND variable from file %s is not 0.\nThe cal pulse condition code will be set to 0 in the GRP variable.\n', ...
            infiles{grp}));
    end
    
  
    GRP.indiv_fnames{grp}=GND.indiv_fnames;
    GRP.indiv_subnames{grp}=GND.indiv_subnames;
    GRP.indiv_traits{grp}=GND.indiv_traits;
    GRP.indiv_art_ics{grp}=GND.indiv_art_ics;
    if ~isempty(GND.cals),
        GRP.cals.indiv_cals{grp}=GND.cals.indiv_cals(use_chans,:,:);
        GRP.cals.indiv_cal_ct{grp}=GND.cals.indiv_cal_ct;
        GRP.cals.grand_cals(:,:,grp)=GND.cals.grand_cals(use_chans,:);
    end
    
    if create_difs,
        bin_sub_ct(:,grp)=GND.sub_ct(use_bins);
        GRP.indiv_bin_ct{grp}=GND.indiv_bin_ct(:,use_bins);
        GRP.indiv_bin_raw_ct{grp}=GND.indiv_bin_raw_ct(:,use_bins);
        % Grab ERPs from whole group and compute sum of squares (SS) for each
        % bin
        group_grands(:,:,:,grp)=GND.grands(use_chans,:,use_bins);
        n_subs=zeros(1,1,n_use_bins);
        n_subs(1,1,:)=GND.sub_ct(use_bins);
        n_subs=repmat(n_subs,[length(use_chans) n_tpts 1]);
        group_ss(:,:,:,grp)=(GND.grands_stder(use_chans,:,use_bins).^2).*n_subs.*(n_subs-1);
    end
end

if create_difs,
    % Compute grand average difference waves & t-scores
    dif_ct=0;
    for grpA=1:n_groups,
        for grpB=(grpA+1):n_groups,
            for b=1:n_use_bins,
                dif_ct=dif_ct+1;
                nA=bin_sub_ct(b,grpA); % # of participants who contributed to this bin from Group A
                nB=bin_sub_ct(b,grpB); % # of participants who contributed to this bin from Group B
                GRP.grands(:,:,dif_ct)=group_grands(:,:,b,grpA)-group_grands(:,:,b,grpB);
                SS_pooled=(group_ss(:,:,b,grpA)+group_ss(:,:,b,grpB))/(nA+nB-2);
                GRP.grands_stder(:,:,dif_ct)=sqrt( SS_pooled*(nA+nB)/(nA*nB) );
                
                % describe difference in bin info
                GRP.bin_info(b).bindesc=sprintf('%s (%s-%s)',bindesc{b}, ...
                    GRP.group_desc{grpA},GRP.group_desc{grpB});
                GRP.bin_info(b).condcode=bin_condcodes(b);
                GRP.bin_info(b).groupA=grpA;
                GRP.bin_info(b).groupB=grpB;
                GRP.bin_info(b).source_binA=b;
                GRP.bin_info(b).source_binB=b;
                GRP.bin_info(b).n_subsA=nA;
                GRP.bin_info(b).n_subsB=nB;
                GRP.bin_info(b).op='A-B';
            end
        end
    end
    GRP.grands_t=GRP.grands./GRP.grands_stder;
end

if isempty(p.Results.out_fname),
    %Create GUI
    [jname, jpath]=uiputfile({'*.GRP','*.GRP files'; ...
        '*.*','All files'},'Save GRP variable as:','untitled.GRP');
    if ~jpath,
        fprintf('Output filename selection cancelled.  GRP variable NOT saved to disk.\n');
    else
        GRP=save_matmk(GRP,jname,jpath,1); % 1 means that user won't be asked again about saving file
    end 
elseif ~strcmpi(p.Results.out_fname,'no save'),
    [jpath, jname]=pathNname(p.Results.out_fname);
    GRP=save_matmk(GRP,jname,jpath);
end
