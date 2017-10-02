%function [pval, t_orig, tmx_ptile, seed_state, est_alpha]=mxt_perm2(dataA,dataB,n_perm,alph,tail,verblevel,seed_state,freq_domain)
%
% mxt_perm2-Two sample permutation test based on a t-statistic and null
% hypothesis of a mean difference of zero.  For one-tailed tests, the null distribution
% is derived from the most extreme t-score of all permutations.  For two-tailed
% tests, the null distribution is derived from the maximum absolute t-score of
% all permutations.  Can handle multiple electrodes and time
% points. According to Hemmelammnn et al. (2004), the t_max permutation
% test has relatively good power for multivariate data whose dimensions are
% highly correlated (like EEG/ERP data).  See Maris (2004) for an
% explanation of how permutation tests like this control for multiple
% comparisons.  Note, Maris uses Hotelling's multivariate T^2 statistic in
% his paper (instead of the t_max used here) and he advocated performing
% the tests over multiple time points (e.g., a 10 ms time window).
% However, I found that the T^2 statistic did not have much power when the
% sample size was near the number of electrodes and the results of the test
% could vary greatly depending on the size of the time window.
%
% Required Inputs:
%  dataA  - 3D matrix of ERPs from Group A (Channel x Time x Participant)
%  dataB  - 3D matrix of ERPs from Group B (Channel x Time x Participant)
%
% Optional Inputs:
%  n_perm - number of permutations {default=1000}.  Manly (1997) suggests
%           using at least 1000 permutations for an alpha level of 0.05 and 
%           at least 5000 permutations for an alpha level of 0.01.
%  alph   - test alpha level {default=.05}
%  tail   - [1 | 0 | -1] If tail=1, the alternative hypothesis is that the
%           means of Group A are greater than those of Group B (upper tailed test).  
%           If tail=0, the alternative hypothesis is that the means of the 
%           two groups are different (two tailed test).  If tail=-1, the 
%           alternative hypothesis is that the means of Group A are less 
%           than those of Group B (lower tailed test). {default: 0}
%  verblevel   - An integer specifiying the amount of information you want
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
% Outputs:
%  pval       - p-value at each time point and electrode (corrected for
%                multiple comparisons via permutation test)
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
% Author:
% David Groppe
% March, 2009
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
% 3/26/2010-fixed critical t-scores (tmx_ptiles) for two-tailed test.  I was
% allowing them to be asymmetrical but this differed from how I was
% computing p-values.  It turns out the way I was computing p-values (from
% Hemmelmann) was effectively using the absolute value of all the t-scores.
% Thus the critical t scores should be symmetrical.
%
% 4/1/2010-Estimated alpha level added
%
% 12/14/2010-freq_domain optional input added

%This code has an appropriate false positive rate when run on simulated
%data at one time point:
%
% tested with script: sim_mxt_perm2.m
% n_sim=4000, n_perm=3000, data=2x3, 16 participants per group, Gaussian, alpha=.05
% Two-tailed test: Empirical alpha rate: 0.0495 (SE=.0034)
% Upper tailed test: Empirical alpha rate: 0.0503 (SE=0.0035)
% Lower tailed test: Empirical alpha rate: 0.0500 (SE=.0034)
%
% n_sim=4000, n_perm=3000, data=2x3, 16 participants Group A, 20 Group B, Gaussian, alpha=.05
% Lower tailed test: Empirical alpha rate: 0.0502 (SE=.0035)
% Upper tailed test: Empirical alpha rate: 0.0483 (SE=.0034)

%%%%%%%%%%%%%%%% FUTURE IMPROVEMENTS %%%%%%%%%%%%%%%%%
%
% -Test for equality of trials or variance?
% -Add alternate test statistics?

function [pval, t_orig, tmx_ptile, seed_state, est_alpha]=mxt_perm2(dataA,dataB,n_perm,alph,tail,verblevel,seed_state,freq_domain)

if nargin<2,
    error('You need to provide data for two groups of participants.');
end

if nargin<3,
    n_perm=2000;
end

if nargin<4,
    alph=.05;
elseif (alph>=1) || (alph<=0)
    error('Argument ''alph'' needs to be between 0 and 1.');
end

if alph<=.01 && n_perm<5000,
    watchit(sprintf('You are probably using too few permutations for an alpha level of %f. Type ">>help mxt_perm2" for more info.',alph));
elseif alph<=.05 && n_perm<1000,
    watchit(sprintf('You are probably using too few permutations for an alpha level of %f. Type ">>help mxt_perm2" for more info.',alph));
end

if nargin<5,
    tail=0;
elseif (tail~=0) && (tail~=1) && (tail~=-1),
    error('Argument ''tail'' needs to be 0,1, or -1.');
end

if nargin<6,
    verblevel=2;
end

if 0,
    defaultStream=RandStream.getDefaultStream; %random # generator state
    if (nargin<7) || isempty(seed_state),
        %Store state of random number generator
        seed_state=defaultStream.State;
    else
        defaultStream.State=seed_state; %reset random number generator to saved state
    end
else
    fprintf('Not seeding random number generator.\n');
    seed_state=NaN;
end

if (nargin<8),
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
    fprintf('mxt_perm2: Number of channels: %d\n',n_chan);
    if freq_domain,
        fprintf('mxt_perm2: Number of frequencies: %d\n',n_tpt);
    else
        fprintf('mxt_perm2: Number of time points: %d\n',n_tpt);
    end
    fprintf('mxt_perm2: Total # of comparisons: %d\n',n_chan*n_tpt);
    fprintf('mxt_perm2: Number of participants in Group A: %d\n',n_subsA);
    fprintf('mxt_perm2: Number of participants in Group B: %d\n',n_subsB);
    fprintf('t-score degrees of freedom: %d\n',total_subs-2); %edit for alt t-test??
end
VerbReport(sprintf('Executing permutation test with %d permutations...',n_perm),2,verblevel);
if (verblevel>=2),
    fprintf('Permutations completed: ');
end

% Factors that are used to compute t-scores.  Saves time to compute them
% now rather than to compute them anew for each permutation.
df=n_subsA+n_subsB-2;
mult_fact=(n_subsA+n_subsB)/(n_subsA*n_subsB);

mx_t=zeros(1,n_perm);
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
    %compute most extreme t-score
    mx_t(perm)=tmax2(data,grp1,grp2,n_subsA,n_subsB,df,mult_fact);
end
    
%End of permutations, print carriage return if it hasn't already been done
%(i.e., perm is NOT a multiple of 1000)
if (verblevel>=2) && rem(perm,1000)
    fprintf('\n');
end

%Compute critical t's
if tail==0,
    %two-tailed, test statistic is biggest absolute value of all t's
    mx_t=abs(mx_t);
    tmx_ptile(2)=prctile(mx_t,100-100*alph);
    tmx_ptile(1)=-tmx_ptile(2);
    est_alpha=mean(mx_t>=tmx_ptile(2));
elseif tail==1,
    %upper tailed
    tmx_ptile=prctile(mx_t,100-100*alph);
    est_alpha=mean(mx_t>=tmx_ptile);
else
    %tail=-1, lower tailed
    tmx_ptile=prctile(mx_t,alph*100);
    est_alpha=mean(mx_t<=tmx_ptile);
end
if verblevel~=0,
    fprintf('Desired family-wise alpha level: %f\n',alph);
    fprintf('Estimated actual family-wise alpha level: %f\n',est_alpha);
end

%Compute t-scores of actual observations
[dummy t_orig]=tmax2(data,1:n_subsA,(n_subsA+1):total_subs,n_subsA,n_subsB,df,mult_fact);
pval=zeros(n_chan,n_tpt);
%compute p-values
for t=1:n_tpt,
    for c=1:n_chan,
        if tail==0,
            pval(c,t)=mean(mx_t>=abs(t_orig(c,t))); %note mx_t are now all positive due to abs command above
        elseif tail==1,
            pval(c,t)=mean(mx_t>=t_orig(c,t));
        elseif tail==-1,
            pval(c,t)=mean(mx_t<=t_orig(c,t));
        end
    end
end


function [mx_t, all_t]=tmax2(dat,grp1,grp2,n_subsA,n_subsB,df,mult_fact)
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
[mx1 id_row]=max(abs(all_t));
[mx2 id_col]=max(mx1);
mx_t=all_t(id_row(id_col),id_col);


