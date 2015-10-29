% rm_bins() - Removes bin from a Mass Univariate ERP Toolbox GND or GRP variable
%             
% Usage:
%  >> GNDorGRP=rm_bins(GNDorGRP,bins)
%
% Required Inputs:
%   GNDorGRP - A  GND or GRP struct variable. GND variables are 
%              produced by the functions avgs2GND.m or sets2GND.m.  GRP 
%              variables are produced by GNDs2GRP.m.  
%   bins     - [integer vector] The bin or bins you would like to remove
%              from the GND or GRP variable.
%
% Output:
%   GNDorGRP - The GND or GRP struct variable with the specified bins and
%              associated permutation test removed.
%
% Example:
%   To remove Bin #41 from a GND variable:
% >> GND=rm_bins(GND,41);
%
%   To remove Bins #3-30 and Bins #40-50 from a GND variable:
% >> GND=rm_bins(GND,[3:30 40:50]);
%
%   To remove Bin #1 from a GRP variable
% >> GRP=rm_bins(GRP,41);
%
% Author:
% David Groppe
% Kutaslab, 4/2010

function GNDorGRP=rm_bins(GNDorGRP,bins)

fldnms=fieldnames(GNDorGRP);
isGRP=0;
if ismember('group_desc',fldnms),
   isGRP=1; 
end

n_bins=length(GNDorGRP.bin_info);
not_have=setdiff(bins,1:n_bins);
if ~isempty(not_have),
    watchit(['This GND or GRP variable does NOT have the following bin(s): ' int2str(not_have)]);
end
bins=intersect(bins,1:n_bins);
if isempty(bins),
    fprintf('Not removing any bins.\n');
else
    fprintf('Removing the following bin(s): %s\n',int2str(bins));
    
    %find any permutation tests performed at that bin and remove them
    n_ptests=length(GNDorGRP.t_tests);
    rm_ptest=[];
    for a=1:n_ptests,
        if ismember(GNDorGRP.t_tests(a).bin,bins)
            rm_ptest=[rm_ptest a];
        end
    end
    if ~isempty(rm_ptest),
        fprintf('Removing the following permutation test results because they were performed on bins being deleted: %s\n',int2str(rm_ptest))
    end
    use_ptests=setdiff(1:n_ptests,rm_ptest);
    GNDorGRP.t_tests=GNDorGRP.t_tests(use_ptests);
    
    use_bins=setdiff(1:n_bins,bins);
    
    GNDorGRP.grands=GNDorGRP.grands(:,:,use_bins);
    GNDorGRP.grands_stder=GNDorGRP.grands_stder(:,:,use_bins);
    GNDorGRP.grands_t=GNDorGRP.grands_t(:,:,use_bins);
    GNDorGRP.bin_info=GNDorGRP.bin_info(use_bins);
    
    if isGRP,
        %fields unique to GRP variables
        n_group=length(GNDorGRP.group_desc);
        for a=1:n_group,
            GNDorGRP.indiv_bin_ct{a}=GNDorGRP.indiv_bin_ct{a}(:,use_bins);
            GNDorGRP.indiv_bin_raw_ct{a}=GNDorGRP.indiv_bin_raw_ct{a}(:,use_bins);
        end
    else
        %fields unique to GND variables
        GNDorGRP.sub_ct=GNDorGRP.sub_ct(use_bins);
        GNDorGRP.indiv_bin_ct=GNDorGRP.indiv_bin_ct(:,use_bins);
        GNDorGRP.indiv_bin_raw_ct=GNDorGRP.indiv_bin_raw_ct(:,use_bins);
        GNDorGRP.indiv_erps=GNDorGRP.indiv_erps(:,:,use_bins,:);
    end
end