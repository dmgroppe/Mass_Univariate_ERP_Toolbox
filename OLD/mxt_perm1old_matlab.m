%function [pval, t_orig, tmx_ptile, seed_state, est_alpha]=mxt_perm1(data,n_perm,alph,tail,verblevel,seed_state,freq_domain)
%
% mxt_perm1-One sample permutation test based on a t-statistic and null
% hypothesis of a mean of zero.  Can handle multiple electrodes and time
% points/frequencies. According to Hemmelammnn et al. (2004), the t_max permutation
% test has relatively good power for multivariate data whose dimensions are
% highly correlated (like EEG/ERP data).  See Maris (2004) for an
% explanation of how permutation tests like this control for multiple
% comparisons.  Note, Maris uses Hotelling's multivariate T^2 statistic in
% his paper (instead of the t_max used here) and he advocated performing
% the tests over multiple time points (e.g., a 10 ms time window).
% However, our lab found that the T^2 statistic did not have much power when the
% sample size was near the number of electrodes and the results of the test
% could vary greatly depending on the size of the time window.
%
% Inputs:
%  data   - 3D matrix of data (Channel x Time x Participant)
%
% Optional Inputs:
%  n_perm      - number of permutations {default=2000}.  Manly (1997) suggests
%                using at least 1000 permutations for an alpha level of 0.05 and 
%                at least 5000 permutations for an alpha level of 0.01.
%  alph        - test alpha level {default=.05}
%  tail        - [1 | 0 | -1] If tail=1, the alternative hypothesis is that the
%                mean of the data is greater than 0 (upper tailed test).  If tail=0,
%                the alternative hypothesis is that the mean of the data is different
%                than 0 (two tailed test).  If tail=-1, the alternative hypothesis
%                is that the mean of the data is less than 0 (lower tailed test).
%                {default: 0}
%  verblevel   - An integer specifying the amount of information you want
%                this function to provide about what it is doing during runtime.
%                 Options are:
%                    0 - quiet, only show errors, warnings, and EEGLAB reports
%                    1 - stuff anyone should probably know
%                    2 - stuff you should know the first time you start working
%                        with a data set {default value}
%                    3 - stuff that might help you debug (show all
%                        reports)
%  seed_state  - The initial state of the random number generating stream
%                (see MATLAB documentation for "randstream"). If you pass
%                a value from a previous run of this function, it should
%                reproduce the exact same values.
%  freq_domain - If 0, command line report will be given in temporal units
%                (e.g. time points).  Otherwise, the report will be given 
%                in frequency domain units (e.g., frequencies). {default:
%                0}
%  
%
% Outputs:
%  pval       - p-value at each time point and electrode (corrected for
%               multiple comparisons via permutation test)
%  t_orig     - t-score at each time point and electrode
%  tmx_ptile  - critical t-scores for given alpha level
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
% Author:
% David Groppe
% June, 2009
% Kutaslab, San Diego
%
% References:
%
% Hemmelmann, et al. (2004) Multivariate tests for the evaluation of
% high-dimensional EEG data. Journal of Neuroscience Methods.
%
% Maris, E. (2004) Randomization tests for ERP topographies and whole 
% spatiotemporal data matrices. Psychophysiology
%
% Manly, B.F.J. (1997) Randomization, bootstrap, and Monte Carlo methods in
% biology. 2nd ed. Chapmn and Hall, London.

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%
% 3/8/2010-rand_state optional input and output added
%
% 3/18/2010-modified to be much faster (time for-loop eliminated, some numbers needed to calculate
% t are computed once for all t-tests; "computational formula" used instead of
% conceptual formula).  New version gives exact same results as old version
%
% 4/1/2010-Estimated alpha level added
%
% 12/10/2010-Differences in the steps used for computing t-scores between
% permutations and observed data could lead to rounding errors that would 
% produce inaccurate p-values for very small sample sizes (e.g., 3 
% participants).  That's been fixed.
%
% 12/14/2010-freq_domain optional input added

%This code has an appropriate false positive rate when run on simulated
%data at one time point (use sim_test1.m):
%
% n_sim=4000, n_perm=2000, data=2x1x16 Gaussian, alpha=.05
% Two-tailed test: Empirical alpha rate: 0.047500
% Upper tailed test: Empirical alpha rate: 0.048500
% Lower tailed test: Empirical alpha rate: 0.045500

function [pval, t_orig, tmx_ptile, seed_state, est_alpha]=mxt_perm1(data,n_perm,alph,tail,verblevel,seed_state,freq_domain)

if nargin<1,
    error('You need to provide data.');
end

if nargin<2,
    n_perm=2000;
end

if nargin<3,
    alph=.05;
elseif (alph>=1) || (alph<=0)
    error('Argument ''alph'' needs to be between 0 and 1.');
end

if alph<=.01 && n_perm<5000,
    watchit(sprintf('You are probably using too few permutations for an alpha level of %f. Type ">>help mxt_perm1" for more info.',alph));
elseif alph<=.05 && n_perm<1000,
    watchit(sprintf('You are probably using too few permutations for an alpha level of %f. Type ">>help mxt_perm1" for more info.',alph));
end

if nargin<4,
    tail=0;
elseif (tail~=0) && (tail~=1) && (tail~=-1),
    error('Argument ''tail'' needs to be 0,1, or -1.');
end

if nargin<5,
    verblevel=2;
end

if 0 % ?? for old versions of Matlab
    defaultStream=RandStream.getDefaultStream; %random # generator state
    if (nargin<6) || isempty(seed_state),
        %Store state of random number generator
        seed_state=defaultStream.State;
    else
        defaultStream.State=seed_state; %reset random number generator to saved state
    end
else
   fprintf('Not seeding random number generator.\n'); 
   seed_state=NaN;
end

if (nargin<7),
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

if verblevel~=0,
    fprintf('mxt_perm1: Number of channels: %d\n',n_chan);
    if freq_domain,
        fprintf('mxt_perm1: Number of frequencies: %d\n',n_pts);
    else
        fprintf('mxt_perm1: Number of time points: %d\n',n_pts);
    end
    fprintf('mxt_perm1: Total # of comparisons: %d\n',n_pts*n_chan);
    fprintf('mxt_perm1: Number of participants: %d\n',n_subs);
    fprintf('t-score degrees of freedom: %d\n',n_subs-1);
end


VerbReport(sprintf('Executing permutation test with %d permutations...',n_perm),2,verblevel);
if (verblevel>=2),
    fprintf('Permutations completed: ');
end


%Constant factor for computing t, speeds up computing t to precalculate
%now
sqrt_nXnM1=sqrt(n_subs*(n_subs-1));

mxt=zeros(1,n_perm*2);
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
    
    %get most extreme t-score (sign isn't immportant since we asumme
    %symmetric distribution of null hypothesis for one sample test)
    mxt(perm)=max(max(abs(t)));
end
mxt(n_perm+1:2*n_perm)=-mxt(1:n_perm); %add the negative of all values since we assumme 
%null hypothesis distribution is symmetric

%End permutations completed line
if (verblevel>=2) && rem(perm,1000)
    fprintf('\n');
end

if tail==0,
    %two-tailed
    tmx_ptile(1)=prctile(mxt,100*alph/2);
    tmx_ptile(2)=-tmx_ptile(1);
    est_alpha=mean(mxt>=tmx_ptile(2))*2;
elseif tail==1,
    %upper tailed
    tmx_ptile=prctile(mxt,100-100*alph);
    est_alpha=mean(mxt>=tmx_ptile);
else
    %tail=-1, lower tailed
    tmx_ptile=prctile(mxt,alph*100);
    est_alpha=mean(mxt<=tmx_ptile);
end
if verblevel~=0,
    fprintf('Desired family-wise alpha level: %f\n',alph);
    fprintf('Estimated actual family-wise alpha level: %f\n',est_alpha);
end

%computes t-scores of observations at all channels and time
%points/frequencies
sm=sum(data,3);
mn=sm/n_subs;
sm_sqrs=sum(data.^2,3)-(sm.^2)/n_subs;
stder=sqrt(sm_sqrs)/sqrt_nXnM1;
t_orig=mn./stder;

%compute p-values
pval=zeros(n_chan,n_pts);
for t=1:n_pts,
    for c=1:n_chan,
        if tail==0,
            pval(c,t)=mean(mxt>=abs(t_orig(c,t)))*2;
        elseif tail==1,
            pval(c,t)=mean(mxt>=t_orig(c,t));
        elseif tail==-1,
            pval(c,t)=mean(mxt<=t_orig(c,t));
        end
    end
end

