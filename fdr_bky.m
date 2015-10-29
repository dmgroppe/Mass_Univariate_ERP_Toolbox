% fdr_bky() - Executes the "two-stage" Benjamini, Krieger, & Yekutieli (2006)
%             procedure for controlling the false discovery rate (FDR) of a
%             family of hypothesis tests. FDR is the expected proportion of
%             rejected hypotheses that are mistakenly rejected (i.e., the null
%             hypothesis is actually true for those tests). FDR is a
%             somewhat less conservative/more powerful method for correcting
%             for multiple comparisons than procedures like Bonferroni
%             correction that provide strong control of the family-wise
%             error rate (i.e., the probability that one or more null
%             hypotheses are mistakenly rejected).
%               The procedure implemented by this function is more powerful
%             than the original Benjamini & Hochberg (1995) procedure when 
%             a considerable percentage of the hypotheses in the family are
%             false. To the best of my knowledge, this procedure is only
%             guaranteed to control FDR if the tests are independent.
%             However, simulations suggest that it can control FDR even
%             when the tests are positively correlated (Benjamini et al., 
%             2006).
%
% Usage:
%  >> [h, crit_p]=fdr_bky(pvals,q,report);
%
% Required Input:
%   pvals - A vector or matrix (two dimensions or more) containing the
%           p-value of each individual test in a family of tests.
%
% Optional Inputs:
%   q       - The desired false discovery rate. {default: 0.05}
%   report  - ['yes' or 'no'] If 'yes', a brief summary of FDR results are
%             output to the MATLAB command line {default: 'no'}
%
%
% Outputs:
%   h       - A binary vector or matrix of the same size as the input "pvals."
%             If the ith element of h is 1, then the test that produced the 
%             ith p-value in pvals is significant (i.e., the null hypothesis
%             of the test is rejected).
%   crit_p  - All p-values less than or equal to crit_p are significant
%             (i.e., their null hypotheses are rejected).  If no p-values are
%             significant, crit_p=0.
%
%
% References:
%   Benjamini, Y., Krieger, A.M., & Yekutieli, D. (2006) Adaptive linear 
%     step-up procedures that control the false discovery rate. Biometrika.
%     93(3), 491-507.
%
%   Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery
%     rate: A practical and powerful approach to multiple testing. Journal
%     of the Royal Statistical Society, Series B (Methodological). 57(1),
%     289-300.
%
% Example:
%   [dummy p_null]=ttest(randn(12,15)); %15 tests where the null hypothesis
%                                       %is true
%   [dummy p_effect]=ttest(randn(12,5)+1); %5 tests where the null
%                                          %hypothesis is false
%   [h crit_p]=fdr_bky([p_null p_effect],.05,'yes');
%
%
% For a review on false discovery rate control and other contemporary
% techniques for correcting for multiple comparisons see:
%
%   Groppe, D.M., Urbach, T.P., & Kutas, M. (2011) Mass univariate analysis 
% of event-related brain potentials/fields I: A critical tutorial review. 
% Psychophysiology, 48(12) pp. 1711-1725, DOI: 10.1111/j.1469-8986.2011.01273.x 
% http://www.cogsci.ucsd.edu/~dgroppe/PUBLICATIONS/mass_uni_preprint1.pdf
%
% Author:
% David M. Groppe
% Kutaslab
% Dept. of Cognitive Science
% University of California, San Diego
% March 25, 2010

function [h crit_p]=fdr_bky(p_values,q,report)

if nargin<1,
    error('You need to provide a vector or matrix of p-values.');
else
    if ~isempty(find(p_values<0,1)),
        error('Some p-values are less than 0.');
    elseif ~isempty(find(p_values>1,1)),
        error('Some p-values are greater than 1.');
    end
end

if nargin<2,
    q=.05;
end

if nargin<3,
    report='no';
end

s=size(p_values);
if (length(s)>2) || s(1)>1,
    p_sorted=sort(reshape(p_values,1,prod(s)));
else
    %p-values are already a row vector
    p_sorted=sort(p_values);
end
m=length(p_sorted); %number of tests

%STEP 1: Run classic Benjamini-Hochberg linear step up FDR procedure (BH) with
%slightly more conservative q value
q_prime=q/(1+q);
[hh crit_p]=fdr_bh(p_sorted,q_prime,'pdep');
r1=sum(hh);

if r1==0,
    %NO hypotheses rejected, stop here
    crit_p=0;
    h=p_values*0;
elseif r1==m,
    %ALL hypotheses rejected, stop here
    crit_p=p_sorted(end); %critical p-value is biggest p-value
    h=p_values<=crit_p;
else
    %Continue on.....
    %STEP 2: r1=Estimated # of false hypotheses
    m0hat=m-r1; % m0hat=estimated # of true hypotheses
    %repeat BH with new q
    q_star=q_prime*m/m0hat;
    [hh crit_p]=fdr_bh(p_sorted,q_star,'pdep');
    h=p_values<=crit_p;
end

if strcmpi(report,'yes'),
    n_sig=sum(hh);
    if n_sig==1,
        fprintf('Out of %d tests, %d is significant using a false discovery rate of %f.\n',m,n_sig,q);
    else
        fprintf('Out of %d tests, %d are significant using a false discovery rate of %f.\n',m,n_sig,q);
    end
    fprintf('FDR two-stage procedure used is guaranteed valid for independent tests.\n');
end


function [h crit_p]=fdr_bh(pvals,q,method)
% fdr_bh() - Executes the Benjamini & Hochberg (1995) procedure for
%            controlling the false discovery rate (FDR) of a family of 
%            hypothesis tests. FDR is the expected proportion of rejected
%            hypotheses that are mistakenly rejected (i.e., the null
%            hypothesis is actually true for those tests). FDR is a
%            somewhat less conservative/more powerful method for correcting
%            for multiple comparisons than procedures like Bonferroni
%            correction that provide strong control of the family-wise
%            error rate (i.e., the probability that one or more null
%            hypotheses are mistakenly rejected).
%
% Usage:
%  >> [h, crit_p]=fdr_bh(pvals,q,method);
%
% Required Input:
%   pvals - A vector or matrix (two dimensions or more) containing the
%           p-value of each individual test in a family of tests.
%   q     - The desired false discovery rate.
%
% Optional Inputs:
%   method  - ['pdep' or 'dep'] If 'pdep,' the original Bejnamini & Hochberg
%             FDR procedure is used, which is guaranteed to be accurate if
%             the individual tests are independent or positively dependent
%             (e.g., positively correlated).  If 'dep,' the FDR procedure
%             described in Benjamini & Yekutieli (2001) that is guaranteed
%             to be accurate for any test dependency structure (e.g.,
%             positively and/or negatively correlated tests) is used. 'dep'
%             is always appropriate to use but is less powerful than 'pdep.'
%             {default: 'pdep'}
%
% Outputs:
%   h       - A binary vector or matrix of the same size as the input "pvals."
%             If the ith element of h is 1, then the test that produced the 
%             ith p-value in pvals is significant (i.e., the null hypothesis
%             of the test is rejected).
%   crit_p  - All p-values less than or equal to crit_p are significant
%             (i.e., their null hypotheses are rejected).  If no p-values are
%             significant, crit_p=0.
%
%
% References:
%   Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery
%     rate: A practical and powerful approach to multiple testing. Journal
%     of the Royal Statistical Society, Series B (Methodological). 57(1),
%     289-300.
%
%   Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery
%     rate in multiple testing under dependency. The Annals of Statistics.
%     29(4), 1165-1188.
%
% Example:
%   [dummy p_null]=ttest(randn(12,15)); %15 tests where the null hypothesis
%                                       %is true
%   [dummy p_effect]=ttest(randn(12,5)+1); %5 tests where the null
%                                          %hypothesis is false
%   [h crit_p]=fdr_bh([p_null p_effect],.05,'pdep','yes');
%
%
% Author:
% David M. Groppe
% Kutaslab
% Dept. of Cognitive Science
% University of California, San Diego
% March 24, 2010


if nargin<3,
    method='pdep';
end

%Note: pvals is already sorted, and a row vector
m=length(pvals); %number of tests

if strcmpi(method,'pdep'),
    %BH procedure for independence or positive dependence
    thresh=[1:m]*q/m;
elseif strcmpi(method,'dep')
    %BH procedure for any dependency structure
    denom=m*sum(1./[1:m]);
    thresh=[1:m]*q/denom;
else
    error('Argument ''method'' needs to be ''pdep'' or ''dep''.');
end

rej=pvals<=thresh;
max_id=find(rej,1,'last'); %find greatest significant pvalue
if isempty(max_id),
    crit_p=0;
    h=pvals*0;
else
    crit_p=pvals(max_id);
    h=pvals<=crit_p;
end




