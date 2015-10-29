function infiles=get_set_infiles(gui_infiles_or_tmplt,sub_ids)
%function infiles=get_set_infiles(gui_infiles_or_tmplt,sub_ids)

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 4/4/2011: Made compatible with Windows.
%

if strcmpi(gui_infiles_or_tmplt,'GUI'),
    loading=1;
    infiles=[];
    while loading,
        [neofname, inpath]=uigetfile({'*.set','*.set files'; ...
            '*clean.set','*clean.set files'; ...
            '*noartIC.set','*noartIC.set files'; ...
            '*.*','All Files (*.*)'},'EEGLAB set Files to Import','MultiSelect','on');
        if ischar(neofname),
            clear infname;
            infname{1}=neofname; %make it a cell array for consistent syntax below
        else
            infname=neofname;
        end
        if ~inpath,
            if isempty(infiles),
                error('File selection cancelled.  Aborting processing.\n');
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
elseif iscell(gui_infiles_or_tmplt)
    infiles=gui_infiles_or_tmplt;
else
    if isempty(sub_ids)
        error('If gui_infiles_or_tmplt is a string filename template, you must specify the participant numbers with the argument ''sub_ids''.');
    elseif ~ismember('#',gui_infiles_or_tmplt),
        error('The filename template %s needs to contain a # to indicate where the participant numbers go (e.g., odbl#.nrm).', ...
            gui_infiles_or_tmplt);
    else
        lb_id=find(gui_infiles_or_tmplt=='#');
        prefix=gui_infiles_or_tmplt(1:(lb_id(1)-1)); %indexing lb_id by 1 in case there are multiple pound signs
        postfix=gui_infiles_or_tmplt((lb_id(1)+1):end);
        n_subs=length(sub_ids);
        infiles=cell(1,n_subs);
        for s=1:n_subs,
            no_pad=[prefix num2str(sub_ids(s)) postfix];
            if sub_ids(s)<10,
                padded=[prefix num2str(sub_ids(s),'%.2d') postfix];
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