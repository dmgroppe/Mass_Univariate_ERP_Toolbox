%function [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm1(data,chan_hood,n_perm,fwer,tail,thresh_p,verblevel,seed_state,freq_domain)
%
% clust_perm1-One sample cluster-based permutation test using the "cluster
%             mass" statistic and a null hypothesis of a mean of zero.  This
%             function can handle multiple electrodes and time points/frequencies.  
%             This test was originally proposed for MRI data by Bullmore et al.
%             (1999) and for EEG/MEG analysis by Maris & Oostenveld (2007).
%
% Required Inputs:
%  data      - 3D matrix of ERP data (Channel x Time x Participant)
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
%  verblevel       - An integer specifying the amount of information you want
%                    this function to provide about what it is doing during runtime.
%                     Options are:
%                      0 - quiet, only show errors, warnings, and EEGLAB reports
%                      1 - stuff anyone should probably know
%                      2 - stuff you should know the first time you start working
%                          with a data set {default value}
%                      3 - stuff that might help you debug (show all
%                        reports)
%  seed_state      - The initial state of the random number generating stream
%                    (see MATLAB documentation for "randstream"). If you pass
%                    a value from a previous run of this function, it should
%                    reproduce the exact same values.
%  freq_domain     - If 0, command line report will be given in temporal units
%                    (e.g. time points).  Otherwise, the report will be given
%                    in frequency domain units (e.g., frequencies). {default:
%                    0}
%
%
% Outputs:
%  pval       - p-value at each time point and electrode (corrected for
%               multiple comparisons via the permutation test)
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
%               impossible to approximately achieve.
%
%
% Notes:
% -Unlike a parametric test (e.g., an ANOVA), a discrete set of p-values
% are possible (at most the number of possible permutations).  Since the
% number of possible permutations grows rapidly with the number of
% participants, this is only an issue for small sample sizes (e.g., 6
% participants).  When you have such a small sample size, the
% limited number of p-values may make the test overly conservative (e.g.,
% you might be forced to use an alpha level of .0286 since it is the biggest
% possible alpha level less than .05).
%
% -For two-tailed tests clust_perm1.m effectively performs an upper-tailed 
% test and a lower-tailed test and multiplies the resulting p-values by 2.
% In other words, it does two tests and then applies Bonferroni correction.
% This is faster than performing a "proper" two-tailed permutation test
% since you only have to form clusters for one polarity (e.g., if you form
% negative clusters for each permutation, the negative of that distribution
% will be valid for computing p-values for postive clusters) but can
% result in p-values greater than 1.
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

%This code has an appropriate false positive rate when run on Gaussian
%noise (26 channels, 25 time points, 16 subs, two-tailed test). Done with
%sim_test_clust1.m
%


function [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm1(data,chan_hood,n_perm,fwer,tail,thresh_p,verblevel,seed_state,freq_domain)

if nargin<1,
    error('You need to provide data.');
end

if nargin<2,
    error('You need to provide a chan_hood matrix.');
end

if nargin<3,
    n_perm=2000;
end

if nargin<4,
    fwer=.05;
elseif (fwer>=1) || (fwer<=0)
    error('Argument ''fwer'' needs to be between 0 and 1.');
end

if fwer<=.01 && n_perm<5000,
    watchit(sprintf('You are probably using too few permutations for a FWER (i.e., alpha level) of %f. Type ">>help clust_perm1" for more info.',fwer));
elseif fwer<=.05 && n_perm<1000,
    watchit(sprintf('You are probably using too few permutations for a FWER (i.e., alpha level) of %f. Type ">>help clust_perm1" for more info.',fwer));
end

if nargin<5,
    tail=0;
elseif (tail~=0) && (tail~=1) && (tail~=-1),
    error('Argument ''tail'' needs to be 0,1, or -1.');
end

if nargin<6,
    thresh_p=.05;
elseif thresh_p<=0 || thresh_p>1,
    error('Argument thresh_p needs to take a value between 0 and 1');
end

if nargin<7,
    verblevel=2;
end

%Get random # generator state
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

if (nargin<9),
    freq_domain=0;
end


s=size(data);
n_chan=s(1);
n_pts=s(2);
n_subs=s(3);
if n_subs<2,
    error('You need data from at least two observations (e.g., participants) to perform a hypothesis test.')
end

if n_subs<7,
    n_psbl_prms=2^n_subs;
    watchit(sprintf(['Due to the very limited number of participants,' ...
        ' the total number of possible permutations is small.\nThus only a limited number of p-values (at most %d) are possible and the test might be overly conservative.'], ...
        n_psbl_prms));
end

df=n_subs-1; %t-score degrees of freedom
if verblevel~=0,
    fprintf('clust_perm1: Number of channels: %d\n',n_chan);
    if freq_domain,
        fprintf('clust_perm1: Number of frequencies: %d\n',n_pts);
    else
        fprintf('clust_perm1: Number of time points: %d\n',n_pts);
    end
    fprintf('clust_perm1: Total # of comparisons: %d\n',n_pts*n_chan);
    fprintf('clust_perm1: Number of participants: %d\n',n_subs);
    fprintf('t-score degrees of freedom: %d\n',df);
end

if tail
    %one tailed test
    thresh_t=tinv(thresh_p,df); %note unless thresh_p is greater than .5, thresh_t will be negative
else
    %two tailed test
    thresh_t=tinv(thresh_p/2,df);
end

VerbReport(sprintf('Executing a cluster-based permutation test with %d permutations...',n_perm),2,verblevel);
if (verblevel>=2),
    fprintf('Permutations completed: ');
end


%Constant factor for computing t, speeds up computing t to precalculate
%now
sqrt_nXnM1=sqrt(n_subs*(n_subs-1));

mn_clust_mass=zeros(1,n_perm);
sn=zeros(1,1,n_subs);
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
    %randomly set sign of each participant's data
    sn(1,1,1:n_subs)=(rand(1,n_subs)>.5)*2-1;
    sn_mtrx=repmat(sn,[n_chan n_pts 1]);
    
    d_perm=data.*sn_mtrx;
    
    %computes t-score of permuted data across all channels and time
    %points or frequencies
    sm=sum(d_perm,3);
    mn=sm/n_subs;
    sm_sqrs=sum(d_perm.^2,3)-(sm.^2)/n_subs;
    stder=sqrt(sm_sqrs)/sqrt_nXnM1;
    t=mn./stder;
    
    [clust_ids, n_clust]=find_clusters(t,thresh_t,chan_hood,-1);

    %get most extremely negative t-score (sign isn't important since we asumme
    %symmetric distribution of null hypothesis for one sample test)
    mn_clust_mass(perm)=find_mn_mass(clust_ids,t,n_clust);
end


%End permutations completed line
if (verblevel>=2) && rem(perm,1000)
    fprintf('\n');
end

% Estimate true FWER of test
if tail==0,
    %two-tailed
    tmx_ptile=prctile(mn_clust_mass,100*fwer/2);
    est_alpha=mean(mn_clust_mass<=tmx_ptile)*2;
else
    %one tailed
    tmx_ptile=prctile(mn_clust_mass,100*fwer);
    est_alpha=mean(mn_clust_mass<=tmx_ptile);
end
if verblevel~=0,
    fprintf('Desired family-wise error rate: %f\n',fwer);
    fprintf('Estimated actual family-wise error rate: %f\n',est_alpha);
end

%computes t-scores of observations at all channels and time
%points/frequencies
sm=sum(data,3);
mn=sm/n_subs;
sm_sqrs=sum(data.^2,3)-(sm.^2)/n_subs;
stder=sqrt(sm_sqrs)/sqrt_nXnM1;
t_orig=mn./stder;


%compute p-values
pval=ones(n_chan,n_pts);
if tail==0,
    %positive clusters
    [clust_ids, n_clust]=find_clusters(t_orig,-thresh_t,chan_hood,1); %note thresh_t is negative by default

    clust_info.pos_clust_pval=ones(1,n_clust);
    clust_info.pos_clust_mass=zeros(1,n_clust);
    clust_info.pos_clust_ids=clust_ids;
    for a=1:n_clust,
        use_ids=find(clust_ids==a);
        clust_mass=sum(t_orig(use_ids)); 
        clust_p=mean(mn_clust_mass<=(-clust_mass))*2; %multiply by 2 since we're effectively doing Bonferroni correcting for doing two tests (an upper tail and lower tail)
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
        clust_p=mean(mn_clust_mass<=clust_mass)*2; %multiply by 2 since we're effectively doing Bonferroni correcting for doing two tests (an upper tail and lower tail)
        pval(use_ids)=clust_p;
        clust_info.neg_clust_pval(a)=clust_p;
        clust_info.neg_clust_mass(a)=clust_mass;
    end
elseif tail==1,
    %upper tailed
    [clust_ids, n_clust]=find_clusters(t_orig,-thresh_t,chan_hood,1); %note thresh_t is negative by default
    clust_info.pos_clust_pval=ones(1,n_clust);
    clust_info.pos_clust_mass=zeros(1,n_clust);
    clust_info.pos_clust_ids=clust_ids;
    for a=1:n_clust,
        use_ids=find(clust_ids==a);
        clust_mass=sum(t_orig(use_ids));
        clust_p=mean(mn_clust_mass<=(-clust_mass));
        pval(use_ids)=clust_p;
        clust_info.pos_clust_pval(a)=clust_p;
        clust_info.pos_clust_mass(a)=clust_mass;
    end
else
    %lower tailed
    [clust_ids, n_clust]=find_clusters(t_orig,thresh_t,chan_hood,-1); %note thresh_t is negative by default
    clust_info.neg_clust_pval=ones(1,n_clust);
    clust_info.neg_clust_mass=zeros(1,n_clust);
    clust_info.neg_clust_ids=clust_ids;
    for a=1:n_clust,
        use_ids=find(clust_ids==a);
        clust_mass=sum(t_orig(use_ids));
        clust_p=mean(mn_clust_mass<=clust_mass);
        pval(use_ids)=clust_p;
        clust_info.neg_clust_pval(a)=clust_p;
        clust_info.neg_clust_mass(a)=clust_mass;
    end
end

%%% End of Main Function %%%


function mn_clust_mass=find_mn_mass(clust_ids,data_t,n_clust)

mn_clust_mass=0;

%looking for most negative cluster mass
for z=1:n_clust,
    use_ids=(clust_ids==z);
    use_mass=sum(data_t(use_ids));
    if use_mass<mn_clust_mass,
        mn_clust_mass=use_mass;
    end
end