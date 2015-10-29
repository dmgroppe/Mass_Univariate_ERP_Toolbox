function erp_var=save_erps(erp_var,new_fname,new_fname_path,force)
% save_erps() - Saves a GND or GRP variable to disk using its current filename 
%               or a newly specified filename.
%             
% Usage:
%  >> erp_var=save_erps(erp_var,new_fname,new_fname_path);
%
% Required Inputs:
%   erp_var - a Mass Univariate ERP Toolbox GND or GRP struct variable
%
% Optional Inputs:
%   new_fname      - ['gui' or a filename]  If the string 'gui,' then a GUI
%                    will be created that you can use to pick a filename.
%                    Otherwise, new_fname should be a string indicating the 
%                    new filname that will be used to save the variable (by
%                    convention this should NOT include the file path and 
%                    should end in the extension "GND" or "GRP"). If not  
%                    specified, the filename and path currently stored in 
%                    the GND/GRP variable will be used.
%   new_fname_path - a string indicating the new filname path that will be 
%                    used to save the variable. If not specified and a new 
%                    filename is specified, the current working directory
%                    will be used.
%   force          - [integer] If 0, user will be prompted to verify
%                    filename before saving. {default: 0}
%
% Outputs:
%   erp_var - Same as input GND/GRP variable, except filename information is
%             updated (if a new filename and/or path were specified) and
%             GND/GRP.saved field is set to 'yes.'  Also, file is saved to 
%             disk.
%
% Note:
% -For some reason, the GUI doesn't work properly on Macs. Specifically, 
% the file selection based on file extension doesn't work.
%
%
% Examples:
% To save a GND variable with the GUI:  
% >>GND=save_erps(GND,'gui');
%
% To save a GND variable to a new filename:
% >>GND=save_erps(GND,'yngvob.GND');
%
% To save a GRP variable using the filename and filepath currently stored 
% in the GRP variable:
% >>GRP=save_erps(GRP);
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
% 4/4/2011: Made compatible with Windows.
%
% 5/24/2011: It turns out the code wasn't fully compatible with Windows. It
% was inserting / instead of \ for Windows file paths.  That problem has
% been fixed.  Thanks to Elisa Filevich for alerting us to the problem (and
% the solution).
%

if nargin<4,
    force=0;
end

try
    %determine what type of variable erp_var is
    fldnames=fieldnames(erp_var);
catch
    error('The first argument of save_erps.m needs to be a GND or GRP variable.');
end
if ismember('group_desc',fldnames),
    var_type='GRP';
else
    var_type='GND';
end

%Path slash convention
if ismac || isunix
    %Mac or Unix OS
    slash='/';
else
    %Windows OS
    slash='\';
end

gui=0;
if nargin>1,
    %More than GND/GRP variable specified.  File name or path will be
    %changed
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
        erp_var.filepath=new_fname_path;
        erp_var.filename=new_fname;
    end
    full_fname=[new_fname_path new_fname]; %not used if GUI called
else
    %Use file name and path already stored in GND/GRP variable 
    if isempty(erp_var.filepath),
        %Path not saved with GND/GRP variable for some reason
        %use current working directory
        erp_var.filepath=[pwd slash];
    end
    new_fname_path=erp_var.filepath;
    new_fname=erp_var.filename;
    full_fname=[new_fname_path new_fname];
end

if ~isempty(full_fname) && ~gui,
    %file name and path specified or taken from GND/GRP variable
    if force,
        resp='y';
    else
        resp=[];
    end
    while ~(strcmpi(resp,'y') || strcmpi(resp,'yes') || strcmpi(resp,'n') || strcmpi(resp,'no'))
        %resp=input(sprintf('Save %s variable as %s? (y or n) ',var_type,full_fname),'s');
        fprintf('Save %s variable as %s? (y or n)',var_type,full_fname);
        resp=input('','s');
        if ~(strcmpi(resp,'y') || strcmpi(resp,'yes') || strcmpi(resp,'n') || strcmpi(resp,'no'))
            fprintf('Please respond with just the letter "n" or the letter "y" (no quotes).\n');
        end
    end
    if strcmpi(resp,'y') || strcmpi(resp,'yes'),
        erp_var.saved='yes';
        if ismac || isunix
            %Mac or Unix OS
            [s ss]=unix(['ls ' new_fname_path]);
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
                GRP=erp_var;
            else
                GND=erp_var;
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
        [jname, jpath]=uiputfile({'*.GRP','*.GRP files'; ...
            '*.*','All files'},sprintf('Save GRP variable as:','untitled.GRP'));
    else
        [jname, jpath]=uiputfile({'*.GND','*.GND files'; ...
            '*.*','All files'},sprintf('Save GND variable as:','untitled.GND'));
    end
    if ~jpath,
        fprintf('Output filename selection cancelled.  %s variable NOT saved to disk.\n',var_type);
    else
        erp_var=save_erps(erp_var,jname,jpath,1); % 1 means that user won't be asked again about saving file
    end
end

