% rm_ttest() - Removes sets of t-test results from a GND or GRP variable
%             
% Usage:
%  >> GNDorGRP=rm_ttest(GNDorGRP,ttest_ids);
%
% Required Inputs:
%   GNDorGRP - A GND or GRP struct variable. GND variables are 
%              produced by the functions avgs2GND.m or sets2GND.m.  GRP 
%              variables are produced by GNDs2GRP.m.  
%   ttest_ids- [integer vector] The index or indices of the set of t-test
%              results you would like to remove from the GND or GRP variable.
%
% Output:
%   GNDorGRP - The GND or GRP struct variable with the specified set of 
%              t-tests removed.
%
% Example:
%   To remove t-Test Set #2 from a GND variable:
% >> GND=rm_ttest(GND,2);
%
%   To remove t-Test Sets #2-4 and #9-11 from a GND variable:
% >> GND=rm_ttest(GND,[2:4 9:11]);
%
%   To remove t-Test Set #1 from a GRP variable
% >> GRP=rm_ttest(GRP,1);
%
% Author:
% David Groppe
% Kutaslab, 4/2010

function GNDorGRP=rm_ttest(GNDorGRP,ttest_ids)

n_ttests=length(GNDorGRP.t_tests);
not_have=setdiff(ttest_ids,1:n_ttests);
if ~isempty(not_have),
    watchit(['This GND or GRP variable does NOT have the following permutation test(s): ' int2str(not_have)]);
end
ttest_ids=intersect(ttest_ids,1:n_ttests);
if isempty(ttest_ids),
    fprintf('Not removing any permutation tests.\n');
else
    fprintf('Removing the following permutation test(s): %s\n',int2str(ttest_ids));
    use_ptests=setdiff(1:n_ttests,ttest_ids);
    GNDorGRP.t_tests=GNDorGRP.t_tests(use_ptests);
end
