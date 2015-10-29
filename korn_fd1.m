% korn_fd1-Korn et al. (2004) permutation based procedure for
%          controlling the "generalized family-wise error rate" (GFWER) of a
%          family of repeated-measure/one-sample t-tests.  This
%          function implements the less computationally intensive/more 
%          conservative GFWER procedure described by Korn and colleagues 
%         (it is called Procedure A* in the appendix of their paper).
%
% Usage:
% >> [p_adj, p_raw, t_raw, seed_state]=korn_fd1(data,n_perm,n_fd,tail,verblevel,seed_state);
%
%
% Required Input:
%  data - A 3D matrix of ERP or difference wave data (channel x time point x
%         participant)
%
% Optional Inputs:
%  n_perm     - Number of permutations to use to estimate the null
%               hypothesis distribution of all possible permutations.
%               {default: 2500}
%  n_fd       - The lower bound on the number of false discoveries
%               whose probability you want to control.  For example, if
%               n_fd is 1 and you reject all null hypotheses with
%               adjusted p-values less than .05, then the probably that the
%               number of false discoveries exceeds 1 is 5% or less.
%               {default: 1}
%  tail       - [1 | 0 | -1] If tail=1, the alternative hypothesis is that the
%               mean of the data is greater than 0 (upper tailed test).  If tail=0,
%               the alternative hypothesis is that the mean of the data is different
%               than 0 (two tailed test).  If tail=-1, the alternative hypothesis
%               is that the mean of the data is less than 0 (lower tailed test).
%               {default: 0}
%  verblevel  - An integer specifying the amount of information you want
%               this function to provide about what it is doing during runtime.
%                 Options are:
%                    0 - quiet, only show errors, warnings, and EEGLAB reports
%                    1 - stuff anyone should probably know
%                    2 - stuff you should know the first time you start working
%                        with a data set {default value}
%                    3 - stuff that might help you debug (show all
%                        reports)
%  seed_state - The initial state of the random number generating stream
%               (see MATLAB documentation for "randstream"). If you pass
%               a value from a previous run of this function, it should
%               reproduce the exact same values.
%
% Outputs:
%  p_adj      - p-values corresponding to each dimension of the hypothesis
%               test after correcting for multiple comparisons
%  p_raw      - p-values corresponding to each dimension of the hypothesis
%               test WITHOUT correction for multiple comparisons
%  t_raw      - t-scores corresponding to each dimension of the hypothesis
%               test
%  seed_state - The initial state of the random number generating stream
%               (see MATLAB documentation for "randstream") used to
%               generate the permutations. You can use this to reproduce
%               the output of this function.
%
% Author:
% David Groppe
% August, 2009
% Kutaslab, San Diego
%
% Reference:
%    Korn, E. L., Troendle, J. F., McShane, L. M., & Simon, R. (2004).
% Controlling the number of false discoveries: Application to high-
% dimensional genomic data. Journal of Statistical Planning and Inference,
% 124, 379-398.
%


%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%
% ?/?/2010-

function [p_adj, p_raw, t_raw, seed_state]=korn_fd1(data,n_perm,n_fd,tail,verblevel,seed_state)

if nargin<1,
    error('You need to provide data.');
end

if nargin<2.
    n_perm=2500;
end

if nargin<3,
    n_fd=1;
end

if nargin<4,
    tail=0;
elseif (tail~=0) && (tail~=1) && (tail~=-1),
    error('Argument ''tail'' needs to be 0,1, or -1.');
end

if nargin<5,
    verblevel=2;
end

if (nargin<6) || isempty(seed_state),
    seed_state=[];
end

[n_chan, n_tpt, n_sub]=size(data);

% if n_sub<7,
%     n_psbl_prms=2^n_sub;
%     watchit(sprintf(['Due to the very limited number of participants,' ...
%         ' the total number of possible permutations is small.\nThus only a limited number of p-values (at most %d) are possible and the test might be overly conservative.'], ...
%     n_psbl_prms));
% end

k=n_chan*n_tpt;  %# of variables

%Step 0 <-Note, Step #'s follow those in Korn et al.'s appendix
count=ones(1,k);

%Step 1
p_tilda=zeros(1,k);

%Perform t-tests on n_perm permutations of the data (Step 2)
[p_raw, t_raw, p_perms, seed_state]=p_perm1(data,n_perm,tail,verblevel,seed_state);

%Sort the p-values from permutations of the data and reshape them to be a
%variable x permutation matrix
p_perms=sort(reshape(p_perms,k,n_perm),1);

%Sort p-values of observations from smallest to largest
[srtd_p, sort_ids]=sort(reshape(p_raw,1,n_chan*n_tpt));
[dummy, unsort_ids]=sort(sort_ids); %indices to return sorted_p to p_raw order

%Pre-allocate memory
p_adj=-ones(1,k);
p_adj_srtd=p_tilda; %should be all zeros

for r=n_fd+1:k,
    %Steps 3-5: Compute p_tilda
    perm_id=n_fd+1;
    count(r)=1+sum(p_perms(perm_id,:)<=srtd_p(r));
    p_tilda(r)=count(r)/(n_perm+1);
        
    %Note: p_tilda(r)=p_adj_srtd(r)=0 for all r<=n_fd
    
    %Step 6: Compute adjusted p-value
    p_adj_srtd(r)=max(p_tilda(n_fd+1:r));
end
p_adj=p_adj_srtd(unsort_ids);
p_adj=reshape(p_adj,n_chan,n_tpt);

