function matmk_var=save_matmk(matmk_var,new_fname,new_fname_path,force)
% save_matmk() - Saves a GND, GRP, specGND, specGRP, or tfrGND variable to  
%                disk using its current filename or a newly specified filename.
%             
% Usage:
%  >> matmk_var=save_matmk(matmk_var,new_fname,new_fname_path);
%
% Required Inputs:
%   matmk_var - a matlabMK GND, GRP, specGND, specGRP, TFR, tfrGND struct 
%               variable
%
% Optional Inputs:
%   new_fname      - ['gui' or a filename]  If the string 'gui,' then a GUI
%                    will be created that you can use to pick a filename.
%                    Otherwise, new_fname should be a string indicating the 
%                    new filname that will be used to save the variable (by
%                    convention this should NOT include the file path and 
%                    should end in the extension "GND", "GRP", "specGND", 
%                    "specGRP", or "tfrGND").  If not specified, the 
%                    filename and path currently stored in the variable
%                    will be used.
%   new_fname_path - a string indicating the new filname path that will be 
%                    used to save the variable. If not specified and a new 
%                    filename is specified, the current working directory
%                    will be used.
%   force          - [integer] If 0, user will be prompted to verify
%                    filename before saving. {default: 0}
%
% Outputs:
%   matmk_var - Same as input GND/GRP/specGND/specGRP/tfrGND variable, except 
%               filename information is updated (if a new filename and/or 
%               path were specified) and GND/GRP/specGND/specGRP.saved 
%               field is set to 'yes.' Also, file is saved to disk.
%
% Note:
% -For some reason the GUI doesn't work properly on Macs. Specifically, 
% the file selection based on file extension doesn't work.
%
%
% Examples:
% To save a GND variable with the GUI:  
% >>GND=save_matmk(GND,'gui');
%
% To save a GND variable to a new filename:
% >>GND=save_matmk(GND,'yngvob.GND');
%
% To save a GRP variable using the filename and filepath currently stored 
% in the GRP variable:
% >>GRP=save_matmk(GRP);
%
% Author:
% David Groppe
% Kutaslab, 3/2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%
% 3/11/10 Added force option, GUI, and updates GND.saved field to 'yes'. 
%
% 3/15/10 Function checks to make sure path exists before attempting save.
%
% 3/23/10 Changed from saveGND.m to save_erps.m as function now works with
% GRP variables too
%
% 12/9/10 Changed from save_erps.m to save_matmk.m as function now works with spectral
% variables too
%
% 1/19/11 Can now handle tfrGND variables
%
% 3/23/11 Can now handle TFR variables 
%
% 5/24/2011 Made compatible with Windows.
%
% 3/5/2012 Gábor Háden made is possible to use path names with spaces in
% the name (when using Linux).


if nargin<4,
    force=0;
end

%Determine what type of variable matmk_var is
var_type=[];
if isfield(matmk_var,'freqs'),
    if isfield(matmk_var,'group_desc'),
        var_type='specGRP';
    else
        var_type='specGND';
    end
elseif isfield(matmk_var,'tfr_fname'),
    var_type='TFR';
elseif isfield(matmk_var,'ftrip'),
    var_type='tfrGND';
elseif isfield(matmk_var,'group_desc'),
    var_type='GRP';
elseif isfield(matmk_var,'grands'),
    var_type='GND';
end
if isempty(var_type),
    error('The first argument of save_matmk.m needs to be a GND, GRP, specGND, specGRP, TFR, or tfrGND variable.');
end

%Figure out which path convention is being used
if ismac || isunix
    %Mac or Unix OS
    slash='/';
else
    %Windows OS
    slash='\';
end

gui=0;
if nargin>1,
    %New filename, GUI, and (possibly) filepath specified
    if strcmpi(new_fname,'gui'),
        gui=1;
    end
    
    %Get Path
    if nargin==2 || isempty(new_fname_path),
        %new filname specified, but not a new path; use pwd
        new_fname_path=[pwd slash];
    end
    if new_fname_path(end)~=slash,
        new_fname_path=[new_fname_path slash];
    end
    if ~gui
        if strcmpi(var_type,'TFR'),
            matmk_var.tfr_fname=[new_fname_path new_fname];
        else
            matmk_var.filepath=new_fname_path;
            matmk_var.filename=new_fname;
        end
    end
    full_fname=[new_fname_path new_fname]; %not used if GUI called
else
    %Use file name and path already stored in matmk variable
    if strcmpi(var_type,'TFR'),
        slash_ids=find(matmk_var.tfr_fname==slash);
        if isempty(slash_ids),
            new_fname_path=[pwd slash];
            full_fname=[new_fname_path matmk_var.tfr_fname];
        else
            new_fname_path=matmk_var.tfr_fname(1:slash_ids(end));
            full_fname=matmk_var.tfr_fname;
        end
    else
        if isempty(matmk_var.filepath),
            %Path not saved with matmk variable for some reason
            %use current working directory
            matmk_var.filepath=[pwd slash];
        end
        new_fname_path=matmk_var.filepath;
        full_fname=[matmk_var.filepath matmk_var.filename];
    end
end

if ~isempty(full_fname) && ~gui,
    %fname specified or taken from GND/GRP/etc... variable
    if force,
        resp='y';
    else
        resp=[];
    end
    while ~(strcmpi(resp,'y') || strcmpi(resp,'yes') || strcmpi(resp,'n') || strcmpi(resp,'no'))
        fprintf('Save %s variable as %s? (y or n)',var_type,full_fname);
        resp=input('','s');
        if ~(strcmpi(resp,'y') || strcmpi(resp,'yes') || strcmpi(resp,'n') || strcmpi(resp,'no'))
            fprintf('Please respond with just the letter "n" or the letter "y" (no quotes).\n');
        end
    end
    if strcmpi(resp,'y') || strcmpi(resp,'yes'),
        matmk_var.saved='yes';
        if ismac || isunix
            %Mac or Unix OS
            [s ss]=unix(['ls "' new_fname_path '"']);
            %Note, s will be 0 if the path exists and 1 otherwise (i.e., 1
            %indicates an error with the command)
        else
            %must be a PC
            if isdir(new_fname_path),
                s=0; %path exists
            else
                s=1; %path doesn't exist
            end
        end
        if s,
            %path doesn't exist
            watchit(sprintf('Path %s does NOT exist.',new_fname_path));
            fprintf('Saving variable to disk aborted. %s variable NOT saved to disk.\n',var_type);
        else
            if strcmpi(var_type,'GRP'),
                GRP=matmk_var;
            elseif strcmpi(var_type,'GND'),
                GND=matmk_var;
            elseif strcmpi(var_type,'specGND'),
                specGND=matmk_var;
            elseif strcmpi(var_type,'specGRP'),
                specGRP=matmk_var;
            elseif strcmpi(var_type,'TFR'),
                TFR=matmk_var;
            else
                tfrGND=matmk_var;
            end
            save(full_fname,var_type);
            fprintf('%s variable successfully saved to disk.\n',var_type);
        end
    else
        fprintf('Saving variable to disk aborted. %s variable NOT saved to disk.\n',var_type);
    end
else
    %Create GUI
    if strcmpi(var_type,'GRP'),
        [jname, jpath]=uiputfile({'*.GRP','*.GRP files',; ...
            '*.*','All files'},sprintf('Save GRP variable as:','untitled.GRP'));
    elseif strcmpi(var_type,'GND'),
        [jname, jpath]=uiputfile({'*.GND','*.GND files'; ...
            '*.*','All files'},sprintf('Save GND variable as:','untitled.GND'));
    elseif strcmpi(var_type,'specGND'),
        [jname, jpath]=uiputfile({'*.specGND','*.specGND files'; ...
            '*.*','All files'},sprintf('Save specGND variable as:','untitled.specGND'));
    elseif strcmpi(var_type,'specGRP'),
        [jname, jpath]=uiputfile({'*.specGRP','*.specGRP files'; ...
            '*.*','All files'},sprintf('Save specGRP variable as:','untitled.specGRP'));
    elseif strcmpi(var_type,'TFR'),
        [jname, jpath]=uiputfile({'*.tfr','*.tfr files'; ...
            '*.*','All files'},sprintf('Save TFR variable as:','untitled.tfr'));
    else
        %tfrGND variable
        [jname, jpath]=uiputfile({'*.tfrGND','*.tfrGND files'; ...
            '*.*','All files'},sprintf('Save tfrGND variable as:','untitled.tfrGND'));
    end
    if ~jpath,
        fprintf('Output filename selection cancelled.  %s variable NOT saved to disk.\n',var_type);
    else
        matmk_var=save_matmk(matmk_var,jname,jpath,1); % 1 means that user won't be asked again about saving file
    end
end

