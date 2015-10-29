% headinfo() - Displays channel, bin and t-test information about the data 
%              stored in a GND or GRP variable (i.e., the Mass Univariate 
%              ERP Toolbox variables for storing ERPs averaged across 
%              participants and between groups respectively) in the MATLAB 
%              command window.  Similar to Kutaslab UNIX "headinfo" command.
%
% Usage:
%  >> headinfo(GND_GRP_or_fname)
%
% Input:
%  GND_GRP_or_fname  - A Mass Univariate ERP Toolbox GND/GRP structure 
%                      variable or the filename of a Mass Univariate ERP
%                      Toolbox GND/GRP structure that has been saved to disk.
%                      To create a GND variable from Kutaslab ERP files (e.g.,
%                      *.nrm files) use avgs2GND.m.  To do the same from
%                      EEGLAB *.set files use sets2GND.m.  To create a
%                      a GRP structure use GNDs2GRP.m. See Mass Univariate ERP 
%                      Toolbox documentation for detailed information about the format
%                      of GND and GRP variables. If you specifiy a filename be
%                      sure to include the file's path, unless the file is
%                      in the current working directory.
%
% Outputs:
%   Text is displayed in command window.
%
% Examples:
% -For GND variable already in memory
% >>headinfo(GND)
%
% -For GND variable on disk but not in MATLAB
% >>headinfo('yngvob.GND')
%
% Author:
% David Groppe
% Kutaslab, 2/2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 3/27/2010-Revised to accommodate GRP variable
%
% 3/30/2010-Revised to deal with GRP/GND variables without bins
%
% 5/6/2010-Revised to deal with FDR control
%
% 6/4/2011-Revised to deal with cluster based tests

function headinfo(GND_GRP_or_fname)

GRP_var=0;
%load GND/GRP variable if a filename was passed
if ischar(GND_GRP_or_fname),
    fprintf('Loading GND or GRP struct from file %s.\n',GND_GRP_or_fname);
    load(GND_GRP_or_fname,'-MAT');
    if ~exist('GND','var') && ~exist('GRP','var')
        error('File %s does not contain a GND or GRP variable.',GND_GRP_or_fname);
    end
    if exist('GRP','var'),
        GND=GRP; %for simplicity GRP variables are re-named GND
        clear GRP;
        GRP_var=1;
    end
else
    GND=GND_GRP_or_fname; %for simplicity GRP and GND variables are named GND
    fldnames=fieldnames(GND);
    if ismember('group_desc',fldnames),
        GRP_var=1;
    end
end

n_chans=length(GND.chanlocs);
n_tpts=length(GND.time_pts);
n_bins=length(GND.bin_info);

fprintf('CHANNELS\n');
for c=1:n_chans,
    fprintf('%d:%s ',c,GND.chanlocs(c).labels);
    if ~mod(c,8)
        fprintf('\n');
    end
end
if mod(c,8),
    fprintf('\n');
end

if GRP_var,
    fprintf('\nGROUPS\n');
    n_groups=length(GND.group_desc);
    for b=1:n_groups,
        fprintf('Group %d: %s (%d participants)\n',b,GND.group_desc{b}, ...
            length(GND.indiv_fnames{b}));
    end
end

fprintf('\nBINS\n');
if n_bins
    for b=1:n_bins,
        fprintf('Bin %d: %s\n',b,GND.bin_info(b).bindesc);
    end
else
    fprintf('No bins are yet stored with these data.\n');
end

%t-Tests
fprintf('\nt-TESTS\n');
n_ptests=length(GND.t_tests);

for t=1:n_ptests,
     fprintf('t-test %d->Bin %d, ',t,GND.t_tests(t).bin);
    if isnan(GND.t_tests(t).n_perm)
        fprintf('q(method=%s)=%g, ',GND.t_tests(t).mult_comp_method, ...
            GND.t_tests(t).desired_alphaORq);
    else
        if strcmpi(GND.t_tests(t).mult_comp_method,'tmax perm test')
            fprintf('Est. Alpha(method=tmax)=%g, ',GND.t_tests(t).estimated_alpha);
        else
            fprintf('Est. Alpha(method=cluster mass)=%g, ',GND.t_tests(t).estimated_alpha);
        end
    end
    if strcmpi(GND.t_tests(t).mean_wind,'yes'),
        fprintf('Mean uV in ');
    else
        fprintf('Each time point in ');
    end
    n_wind=size(GND.t_tests(t).time_wind,1);
    for w=1:n_wind,
        if w==n_wind
            fprintf('%d-%d ms, ',GND.t_tests(t).time_wind(w,1), ...
                GND.t_tests(t).time_wind(w,2));
        else
            fprintf('%d-%d, ',GND.t_tests(t).time_wind(w,1), ...
                GND.t_tests(t).time_wind(w,2));
        end
    end
    ex_ids=setdiff(1:length(GND.chanlocs),GND.t_tests(t).used_chan_ids);
    if isempty(ex_ids),
        fprintf('All channels included\n');
    else
        fprintf('Excluding channels:');
        for c=ex_ids,
            fprintf(' %s',GND.chanlocs(c).labels);
        end
        fprintf('\n');
    end
end
if ~n_ptests,
    fprintf('No t-test results have been stored with these data.\n');
end
