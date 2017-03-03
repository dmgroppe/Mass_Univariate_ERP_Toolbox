%[pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm2(dataA,dataB,chan_hood,n_perm,fwer,tail,thresh_p,verblevel,seed_state,freq_domain)
%
% clust_perm2-Independent samples cluster-based permutation test using the "cluster
%             mass" statistic and a null hypothesis of a mean of zero.  This
%             function can handle multiple electrodes and time points/frequencies.  
%             This test was originally proposed for MRI data by Bullmore et al.
%             (1999) and for EEG/MEG analysis by Maris & Oostenveld (2007).
%
% Required Inputs:
%  dataA      - 3D matrix of ERP data (Channel x Time x Participant)
%  dataB      - 3D matrix of ERP data (Channel x Time x Participant)
%  chan_hood - 2D symmetric binary matrix that indicates which channels are
%              considered neighbors of other channels. E.g., if
%              chan_hood(2,10)=1, then Channel 2 and Channel 10 are
%              nieghbors. You can produce a chan_hood matrix using the
%              function spatial_neighbors.m.
%
% Optional Inputs:
%  n_perm          - Number of permutations {default=2000}.  Manly (1997) suggests
%                    using at least 1000 permutations for an alpha level of 0.05 and
%                    at least 5000 permutations for an alpha level of 0.01.
%  fwer            - Desired family-wise error rate (i.e., alpha level) {default=.05}
%  tail            - [1 | 0 | -1] If tail=1, the alternative hypothesis is that the
%                    mean of the data is greater than 0 (upper tailed test).  If tail=0,
%                    the alternative hypothesis is that the mean of the data is different
%                    than 0 (two tailed test).  If tail=-1, the alternative hypothesis
%                    is that the mean of the data is less than 0 (lower tailed test).
%                    {default: 0}
%  thresh_p        - The test-wise p-value threshold for cluster inclusion. If
%                    a channel/time-point has a t-score that corresponds to an
%                    uncorrected p-value greater than thresh_p, it is assigned
%                    a p-value of 1 and not considered for clustering. Note
%                    that thresh_p automatically takes into account the tail of
%                    the test (e.g., you will get positive and negative t-score
%                    thresholds for a two-tailed test).
%  verblevel       - An integer specifiying the amount of information you want
%                    this function to provide about what it is doing during runtime.
%                     Options are:
%                      0 - quiet, only show errors, warnings, and EEGLAB reports
%                      1 - stuff anyone should probably know
%                      2 - stuff you should know the first time you start working
%                          with a data set {default value}
%                      3 - stuff that might help you debug (show all
%                          reports)
%  seed_state      - The initial state of the random number generating stream
%                    (see MATLAB documentation for "randstream"). If you pass
%                    a value from a previous run of this function, it should
%                    reproduce the exact same values.
%  freq_domain     - If 0, command line report will be given in temporal units
%                    (e.g. time points).  Otherwise, the report will be given 
%                    in frequency domain units (e.g., frequencies). 
%                    {default: 0}
%
% Outputs:
%  pval       - p-value at each time point and electrode (corrected for
%                multiple comparisons via the permutation test)
%  t_orig     - t-score at each time point and electrode
%  clust_info - A struct variable containing information about the
%               clusters found.  Depending on the tail of the test it will
%               be composed of all or half of the following fields:
%                 pos_clust_pval: p-values of the positive clusters
%                 pos_clust_mass: t-score mass of the positive clusters
%                 pos_clust_ids:  channel x time point matrix that
%                   indicated which positive cluster each channel/time point
%                   pair belongs to. E.g., if pos_clust_ids(1,2)=4, then
%                   the second time point of the first channel is a
%                   member of the fourth cluster. 0 indicates that the
%                   channel/time point is not a member of any positive
%                   cluster.
%                 neg_clust_pval: p-values of the negative clusters
%                 neg_clust_mass: t-score mass of the negative clusters
%                 neg_clust_ids:  channel x time point matrix that
%                   indicated which negative cluster each channel/time point
%                   pair belongs to. E.g., if neg_clust_ids(1,2)=4, then
%                   the second time point of the first channel is a
%                   member of the fourth cluster. 0 indicates that the
%                   channel/time point is not a member of any negative
%                   cluster.
%  seed_state - The initial state of the random number generating stream
%               (see MATLAB documentation for "randstream") used to 
%               generate the permutations. You can use this to reproduce
%               the output of this function.
%  est_alpha  - The estimated family-wise alpha level of the test.  With 
%               permutation tests, a finite number of p-values are possible.
%               This function tries to use an alpha level that is as close 
%               as possible to the desired alpha level.  However, if the 
%               sample size is small, a very limited number of p-values are 
%               possible and the desired family-wise alpha level may be 
%               impossible to acheive.
%
% Notes:
% -The null hypothesis of this test is that the data in the two groups are
% "exchangeable" (i.e., the data from a participant in Group A was just as
% likely to have come from Group B).  If the distributions of the two
% populations you're comparing differ in *any* way (i.e., not just in
% central tendency), the null hypothesis will be violated.  For example, if
% the ERPs of the two groups differ in their variability, but not their
% average value (e.g., one set of ERPs is from a group of children and the
% other adults), you'll be more likely than alpha to reject the null
% hypothesis.
%
% -Unlike a parametric test (e.g., an ANOVA), a discrete set of p-values
% are possible (at most the number of possible permutations).  Since the
% number of possible permutations grows rapdily with the number of
% participants, this is only an issue for small sample sizes (e.g., 3
% participants in each group).  When you have such a small sample size, the
% limited number of p-values may make the test overly conservative (e.g., 
% you might be forced to use an alpha level of .0286 since it is the biggest
% possible alpha level less than .05).
%
% -When doing two-tailed tests it is possible to get p-values greater than 1.
% These can be treated equivalent to a p-value of 1. The reason for
% this is that the p-values for positive and negative clusters are computed
% from two different distributions of permuted values.  For a non-cluster
% based permutation test, there is only one distribution of permuted
% values.
%
% Author:
% David Groppe
% May, 2011
% Kutaslab, San Diego
%
% References:
% Bullmore, E. T., Suckling, J., Overmeyer, S., Rabe-Hesketh, S., Taylor, 
% E., & Brammer, M. J. (1999). Global, voxel, and cluster tests, by theory 
% and permutation, for a difference between two groups of structural MR 
% images of the brain. IEEE Transactions on Medical Imaging, 18(1), 32-42. 
% doi:10.1109/42.750253
%
% Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of 
% EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 177-190. 
% doi:10.1016/j.jneumeth.2007.03.024
%
% Manly, B.F.J. (1997) Randomization, bootstrap, and Monte Carlo methods in
% biology. 2nd ed. Chapmn and Hall, London.

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%
% 12/11/2011-Now uses Amy Guthormsen's recursiveless find_clusters.m.
%
% 3/17/2013-Randomization is now compatible with Matlab v13. Thanks to
% Aaron Newman for the fix.
%

%
%This code has an appropriate false positive rate when run on Gaussian
%noise (26 channels, 25 time points, 16 subs, two-tailed test). Done with
%sim_test_clust2.m
%

%%%%%%%%%%%%%%%% FUTURE IMPROVEMENTS %%%%%%%%%%%%%%%%%
%
% -Test for equality of trials or variance?
% -Add alternate test statistics?

function [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm2(dataA,dataB,chan_hood,n_perm,fwer,tail,thresh_p,verblevel,seed_state,freq_domain)

if nargin<2,
    error('You need to provide data for two groups of participants.');
end

if nargin<3,
    error('You need to provide a chan_hood matrix.');
end

if nargin<4,
    n_perm=2000;
end

if nargin<5,
    fwer=.05;
elseif (fwer>=1) || (fwer<=0)
    error('Argument ''fwer'' needs to be between 0 and 1.');
end

if fwer<=.01 && n_perm<5000,
    watchit(sprintf('You are probably using too few permutations for an alpha level of %f. Type ">>help clust_perm2" for more info.',fwer));
elseif fwer<=.05 && n_perm<1000,
    watchit(sprintf('You are probably using too few permutations for an alpha level of %f. Type ">>help clust_perm2" for more info.',fwer));
end

if nargin<6,
    tail=0;
elseif (tail~=0) && (tail~=1) && (tail~=-1),
    error('Argument ''tail'' needs to be 0,1, or -1.');
end

if nargin<7,
    thresh_p=.05;
elseif thresh_p<=0 || thresh_p>1,
    error('Argument thresh_p needs to take a value between 0 and 1');
end

if nargin<8,
    verblevel=2;
end

%get random # generator state
if verLessThan('matlab','7.6')
    watchit('Your version of MATLAB is too old to seed random number generator. You will not be able to exactly reproduce test results.');
    seed_state=NaN;
else
    if verLessThan('matlab','8.1')
        defaultStream=RandStream.getDefaultStream;
    else
        defaultStream=RandStream.getGlobalStream;
    end
    if (nargin<6) || isempty(seed_state),
        %Store state of random number generator
        seed_state=defaultStream.State;
    else
        defaultStream.State=seed_state; %reset random number generator to saved state
    end
end

if (nargin<10),
   freq_domain=0; 
end

if length(size(dataA))~=3 || length(size(dataB))~=3 
    error('dataA and dataB need to be three dimensional (chan x time x participant)')
end
[n_chan, n_tpt, n_subsA]=size(dataA);
[n_chan2, n_tpt2, n_subsB]=size(dataB);

if verblevel,
    warning('off','all'); %for large # of subjects, nchoosek warns that its result is approximate
    n_psbl_prms=nchoosek(n_subsA+n_subsB,n_subsA);
    if n_psbl_prms<100,
        watchit(sprintf(['Due to the very limited number of participants in each group,' ...
            ' the total number of possible permutations is small.\nThus only a limited number of p-values (at most %d) are possible and the test might be overly conservative.'], ...
            n_psbl_prms));
    end
    warning('on','all');
end

if n_chan~=n_chan2,
    error('The number of channels in Group 1 (dataA) and Group 2 (dataB) need to be equal.');
elseif n_tpt~=n_tpt2,
    error('The number of time points in Group 1 (dataA) and Group 2 (dataB) need to be equal.');
end
%combine data
total_subs=n_subsA+n_subsB;
data=dataA; 
data(:,:,(n_subsA+1):total_subs)=dataB;

if verblevel~=0,
    fprintf('clust_perm2: Number of channels: %d\n',n_chan);
    if freq_domain,
        fprintf('clust_perm2: Number of frequencies: %d\n',n_tpt);
    else
        fprintf('clust_perm2: Number of time points: %d\n',n_tpt);
    end
    fprintf('clust_perm2: Total # of comparisons: %d\n',n_chan*n_tpt);
    fprintf('clust_perm2: Number of participants in Group A: %d\n',n_subsA);
    fprintf('clust_perm2: Number of participants in Group B: %d\n',n_subsB);
    fprintf('t-score degrees of freedom: %d\n',total_subs-2);
end
VerbReport(sprintf('Executing permutation test with %d permutations...',n_perm),2,verblevel);
if (verblevel>=2),
    fprintf('Permutations completed: ');
end

% Factors that are used to compute t-scores.  Saves time to compute them
% now rather than to compute them anew for each permutation.
df=n_subsA+n_subsB-2;
mult_fact=(n_subsA+n_subsB)/(n_subsA*n_subsB);
if tail
    %one tailed test
    thresh_t=tinv(thresh_p,df); %note unless thresh_p is greater than .5, thresh_t will be negative
else
    %two tailed test
    thresh_t=tinv(thresh_p/2,df);
end

if tail==0
    mx_clust_massNEG=zeros(1,n_perm);
    mx_clust_massPOS=zeros(1,n_perm);
else
    mx_clust_mass=zeros(1,n_perm);
end
for perm=1:n_perm
    if ~rem(perm,100)
        if (verblevel>=2),
            if ~rem(perm-100,1000)
                fprintf('%d',perm);
            else
                fprintf(', %d',perm);
            end
            if ~rem(perm,1000)
                fprintf('\n');
            end
        end
    end
    %randomly assign participants to conditions
    r=randperm(total_subs);
    grp1=r(1:n_subsA);
    grp2=r((n_subsA+1):total_subs);

    %compute t-scores
    all_t=tmax2(data,grp1,grp2,n_subsA,n_subsB,df,mult_fact);
    
    %form t-scores into clusters
    if tail==0,
        %two-tailed test
        
        %positive clusters
        [clust_ids, n_clust]=find_clusters(all_t,-thresh_t,chan_hood,1); %note, thresh_t should be negative by default
        mx_clust_massPOS(perm)=find_mx_mass(clust_ids,all_t,n_clust,1);
        
        %negative clusters
        [clust_ids, n_clust]=find_clusters(all_t,thresh_t,chan_hood,-1);
        mx_clust_massNEG(perm)=find_mx_mass(clust_ids,all_t,n_clust,-1);
        
    elseif tail>0
        %upper tailed test
        [clust_ids, n_clust]=find_clusters(all_t,-thresh_t,chan_hood,1); %note, thresh_t should be negative by default
        mx_clust_mass(perm)=find_mx_mass(clust_ids,all_t,n_clust,1);
    else
        %lower tailed test
        [clust_ids, n_clust]=find_clusters(all_t,thresh_t,chan_hood,-1); %note, thresh_t should be negative by default
        mx_clust_mass(perm)=find_mx_mass(clust_ids,all_t,n_clust,-1);
    end
end

%End of permutations, print carriage return if it hasn't already been done
%(i.e., perm is NOT a multiple of 1000)
if (verblevel>=2) && rem(perm,1000)
    fprintf('\n');
end

%Compute critical t's
if tail==0,    
    %two-tailed, test statistic is biggest absolute value of all t's
    mx_clust_mass=abs([mx_clust_massPOS mx_clust_massNEG]);
    tmx_ptile(2)=prctile(mx_clust_mass,100-100*fwer);
    tmx_ptile(1)=-tmx_ptile(2);
    est_alpha=mean(mx_clust_mass>=tmx_ptile(2));
elseif tail==1,
    %upper tailed
    tmx_ptile=prctile(mx_clust_mass,100-100*fwer);
    est_alpha=mean(mx_clust_mass>=tmx_ptile);
else
    %tail=-1, lower tailed
    tmx_ptile=prctile(mx_clust_mass,fwer*100);
    est_alpha=mean(mx_clust_mass<=tmx_ptile);
end
if verblevel~=0,
    fprintf('Desired family-wise alpha level: %f\n',fwer);
    fprintf('Estimated actual family-wise alpha level: %f\n',est_alpha);
end

%Compute t-scores of actual observations
t_orig=tmax2(data,1:n_subsA,(n_subsA+1):total_subs,n_subsA,n_subsB,df,mult_fact);

%compute p-values
pval=ones(n_chan,n_tpt);
if tail==0,
    %two-tailed test
    pval=pval*2; %default p-value for channel/time point pairs is 2
    
    %positive clusters
    [clust_ids, n_clust]=find_clusters(t_orig,-thresh_t,chan_hood,1); %note thresh_t is negative by default
    clust_info.pos_clust_pval=ones(1,n_clust);
    clust_info.pos_clust_mass=zeros(1,n_clust);
    clust_info.pos_clust_ids=clust_ids;
    for a=1:n_clust,
        use_ids=find(clust_ids==a);
        clust_mass=sum(t_orig(use_ids));
        clust_p=mean(mx_clust_massPOS>=clust_mass)*2; %multiply by two to correct for doing an upper AND lower tailed test
        pval(use_ids)=clust_p;
        clust_info.pos_clust_pval(a)=clust_p;
        clust_info.pos_clust_mass(a)=clust_mass;
    end
    
    %negative clusters
    [clust_ids, n_clust]=find_clusters(t_orig,thresh_t,chan_hood,-1); %note thresh_t is negative by default
    clust_info.neg_clust_pval=ones(1,n_clust);
    clust_info.neg_clust_mass=zeros(1,n_clust);
    clust_info.neg_clust_ids=clust_ids;
    for a=1:n_clust,
        use_ids=find(clust_ids==a);
        clust_mass=sum(t_orig(use_ids));
        clust_p=mean(mx_clust_massNEG<=clust_mass)*2; %multiply by two to correct for doing an upper AND lower tailed test
        pval(use_ids)=clust_p;
        clust_info.neg_clust_pval(a)=clust_p;
        clust_info.neg_clust_mass(a)=clust_mass;
    end
elseif tail>0
    %positive clusters
    [clust_ids, n_clust]=find_clusters(t_orig,-thresh_t,chan_hood,1); %note thresh_t is negative by default
    clust_info.pos_clust_pval=ones(1,n_clust);
    clust_info.pos_clust_mass=zeros(1,n_clust);
    clust_info.pos_clust_ids=clust_ids;
    for a=1:n_clust,
        use_ids=find(clust_ids==a);
        clust_mass=sum(t_orig(use_ids));
        clust_p=mean(mx_clust_mass>=clust_mass); 
        pval(use_ids)=clust_p;
        clust_info.pos_clust_pval(a)=clust_p;
        clust_info.pos_clust_mass(a)=clust_mass;
    end
else
    %negative clusters
    [clust_ids, n_clust]=find_clusters(t_orig,thresh_t,chan_hood,-1); %note thresh_t is negative by default
    clust_info.neg_clust_pval=ones(1,n_clust);
    clust_info.neg_clust_mass=zeros(1,n_clust);
    clust_info.neg_clust_ids=clust_ids;
    for a=1:n_clust,
        use_ids=find(clust_ids==a);
        clust_mass=sum(t_orig(use_ids));
        clust_p=mean(mx_clust_mass<=clust_mass);
        pval(use_ids)=clust_p;
        clust_info.neg_clust_pval(a)=clust_p;
        clust_info.neg_clust_mass(a)=clust_mass;
    end
end


%%% End of Main Function %%%

function all_t=tmax2(dat,grp1,grp2,n_subsA,n_subsB,df,mult_fact)
% might make this faster by moving it into the code

x1=dat(:,:,grp1);
x2=dat(:,:,grp2);

sm1=sum(x1,3);
mn1=sm1/n_subsA;
ss1=sum(x1.^2,3)-(sm1.^2)/n_subsA;

sm2=sum(x2,3);
mn2=sm2/n_subsB;
ss2=sum(x2.^2,3)-(sm2.^2)/n_subsB;

pooled_var=(ss1+ss2)/df;
stder=sqrt(pooled_var*mult_fact);

all_t=(mn1-mn2)./stder;


function mx_clust_mass=find_mx_mass(clust_ids,data_t,n_clust,tail)

mx_clust_mass=0;
if tail<0
    %looking for most negative cluster mass
    for z=1:n_clust,
        use_ids=(clust_ids==z);
        use_mass=sum(data_t(use_ids));
        if use_mass<mx_clust_mass,
            mx_clust_mass=use_mass;
        end
    end
elseif tail>0,
    %looking for most positive cluster mass
    for z=1:n_clust,
        use_ids=(clust_ids==z);
        use_mass=sum(data_t(use_ids));
        if use_mass>mx_clust_mass,
            mx_clust_mass=use_mass;
        end
    end
end



