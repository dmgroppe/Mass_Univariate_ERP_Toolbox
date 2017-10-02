# Mass_Univariate_ERP_Toolbox

### MATLAB functions for analyzing and visualizing large numbers of *t*-tests performed on event-related potential data.

The Mass Univariate ERP Toolbox is a freely available set of MATLAB functions for performing mass univariate analyses of event-related potentials (ERPs), a noninvasive measure of neural activity popular in cognitive neuroscience. A mass univariate analysis is the analysis of a massive number of simultaneously measured dependent variables via the performance of univariate hypothesis tests (e.g., *t*-tests).  Savvy corrections for multiple comparisons are applied to make spurious findings unlikely while still retaining a useful degree of statistical power. The advantages of mass univariate analyses include:
  * Reduced need for a priori defined time windows/regions of interest
  * Discovery of unexpected effects even when a priori time windows/regions of interest are available
  * Greater spatial and temporal resolution than conventional mean time window analyses

The disadvantages of mass univariate analyses are that they lose some statistical power due to correction for multiple comparisons and some popular corrections for multiple comparisons are not guaranteed to work or may not provide the degree of certainty provided by selective analyses of a priori time windows/regions of interest.  Currently the toolbox supports within-subject and between-subject *t*-tests with false discovery rate controls and control of the family-wise error rate via permutation tests.

This toolbox was produced by members of the [Kutaslab](http://kutaslab.ucsd.edu/) of the [Department of Cognitive Science](http://www.cogsci.ucsd.edu/) at the University of California, San Diego.  If you use the toolbox to perform analyses or to produce figures used in a publication, please cite the following article:

[Groppe, D.M., Urbach, T.P., Kutas, M. (2011) Mass univariate analysis of event-related brain potentials/fields I: A critical tutorial review, *Psychophysiology*, 48(12) pp. 1711-1725, DOI: 10.1111/j.1469-8986.2011.01273.x](http://www.cogsci.ucsd.edu/~dgroppe/PUBLICATIONS/mass_uni_preprint1.pdf)



---
### Documentation and a tutorial for using the code are available here:
[http://openwetware.org/wiki/Mass_Univariate_ERP_Toolbox](http://openwetware.org/wiki/Mass_Univariate_ERP_Toolbox)

---
### Extensions:
Eric Fields has extended this toolbox via his [Factorial Mass Univariate ERP Toolbox](https://github.com/ericcfields/FMUT) (FMUT). The FMUT implements mass univariate statistics for ANOVA-based ERP analyses.

---
DISCLAIMER: The Mass Univariate ERP toolbox is written and released for research purposes only with no guarantee of suitability for any particular purpose. This software, or data obtained from this software, should not under any circumstances be used for clinical purposes.
