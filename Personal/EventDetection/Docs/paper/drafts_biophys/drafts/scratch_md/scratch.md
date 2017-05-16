# Title Page

# Abstract

# Main Text

## Introduction {#s:Intro}

In single-molecule force spectroscopy (SMFS) experiments, a force probe attaches to a molecule and stretches it while measuring force and extension over time (\fRef{Cartoon}). These data are transformed into information such as kinetic rates of processive enzymes,@comstock_direct_2015 protein-ligand bond strength, @yuan_energy_2000 and the energy landscapes of proteins and nucleic acids @dudko_Theory_2008. The location and properties of the *ruptures* in the data (see Figure {#f:Rupture}) are required for many common  analyses, such as applying polymer models and determining the molecular energy landscape.

Atomic force microscopy (AFM) is a powerful tool for studying the mechanical 
properties of molecules.  AFM as an imaging technique can resolve sub-nanometer
 molecular structure such as the major and minor grooves of 
DNA,@ido_beyond_2013 lattice structure of membrane-bound proteins,
@muller_surface_1999 and real-time motion of motor proteins
@ando_high-speed_2013. As a force spectroscopy technique, AFM is capable of
 dynamic experiments such as measuring the unfolding and refolding kinetics 
of single proteins,@he_direct_2015 unzipping double-stranded 
DNA, @krautbauer_Unzipping_2003 and determining the unfolding and refolding 
pathways for membrane proteins @yu_hidden_2017. The viability of AFM in a wide
 range of temperatures, solvents, and other environmental variables makes it 
attractive for studying biological systems. 

During an SMFS-AFM experiment, a cantilever-bound tip with a nanometer-sized radius interacts with a sample (Figure {#f:Cartoon}). The interaction is measured via an optical lever arm system @meyer_novel_1998 in which the displacement of the tip is recorded via deflection of a laser focused on the cantilever. A calibrated tip measures interaction forces from piconewtons to nanonewtons. 

![**Figure {#f:flowchart}.** Algorithm flowchart So that
blah blah blah
linebreak
blah blah](figures/flowchart.png)

The analysis and interpretation of SMFS experiments is dependent on the attachment of the probe to the molecule of interest. In AFM, quality attachments between the tip and the sample may occur less than one curve in tens of thousands.@bosshart_reference-free_2012 Strategies exist to improve attachment rates by coating the AFM tip with a molecule that binds readily to the sample of interest.@walder_robert_rapid_nodate However, the majority of data is still uninterpretable and must be removed. Although a trained human is capable of sorting SMFS data and locating possible events (see Figure {#f:Rupture}), this process is time-consuming and not scientifically reproducible. Excluding data without events and identifying the locations of events is a major challenge in SMFS data management and analysis.

Methods exist to automate detecting events within SMFS data. Automation removes the burden of manual data filtering and improves scientific reproducibility. Techniques for automated event detection in \fec{}s include aligning by the contour length at each extension by applying polymer models;@bosshart_reference-free_2012,kuhn_automated_2005 thresholding based on signal or noise characteristics;@gergely_semi-automatized_2001,roduit_openfovea:_2012 and classification based on transformations of the data into frequency or derivative spaces @kasas_fuzzy_2000,garcia-masso_automated_2016,benitez_searching_2017. These methods do provide an increased degree of autonomy, but their use is limited by their lack of generalization. Specifically, contour-length alignment algorithms bias results towards dominant features and necessarily require a polymer model for the contour length (as a function of force and extension). In SMFS studies of molecules with no existing polymer model or which exhibit rare behavior, alignment algorithms have limited use.  Thresholding and classification-based algorithms generally require optimizing many hard-to-interpret parameters. As shown below, these methods do not generalize well to the typical range of SMFS data (Figure {#f:Performance}).

This thesis describes a new method for detecting events in \fec{}s.  The algorithm, named FEATHER (\acronym{}), requires no \textit{a priori} knowledge of the polymer under study, does not bias data interpretation towards the dominant behavior of the data, and has two easy-to-interpret parameters which generalize well. FEATHER is designed to make few assumptions about the data, operate over a wide range of SMFS experimental conditions, and require a small time footprint compared to existing techniques.

 Chapter 2 of this thesis describes the sample preparation and data acquisition used to measure the \fec{}s of functionalized double-stranded DNA pulled by functionalized AFM cantilevers (Figure {#f:Cartoon}).  The unbinding events in the force-extension curves are then manually annotated at multiple pulling velocities, effective contour lengths, and events per curve (see \tRef{statistics}). Finally, the details and improved performance of FEATHER are described. 


## Materials and Methods

### Sample Preparation

Site-specific chemistry is used to improve the acquisition rate and quality of data. The procedure for surface, sample, and cantilever preparation is described briefly in Appendix A and in detail elsewhere.@walder_robert_rapid_nodate Briefly, through polymerase chain reaction, double-stranded DNA is functionalized with dibenzocyclooctyl (DBCO) at the 5' end of one DNA strand to ensure a covalent bond with an azide-functionalized surface. The DNA is also functionalized with biotin at the other 5' end to ensure a specific but reversible bond with a streptavidin-coated cantilever. These two bonds ensure the tip-DNA bond is broken before the surface-DNA bond, preventing tip contamination. 

### Atomic force microscopy

All atomic force microscopy (AFM) measurements were carried out using an Asylum Cypher Asylum Research}{Cypher ES}. The spring constant and sensitivity were determined using the equipartition theorem method.@hutter_calibration_1993 All measurements were carried out with a 2 second pause at the surface. To promote molecular attachments to the tip, measurements were taken over a series of square, 25}{\textmu{}m} grids at points separated by 1}{\textmu{}m}. The stage was translated to a new grid point after measuring three force-extension curves. 

### Data Annotation

![**Figure {#f:flowchart}.** Algorithm flowchart So that
blah blah blah
linebreak
blah blah](figures/flowchart.png)


Two hundred force-extension curves with events were obtained at three pulling velocities (100 nm/s, 500 nm/s, 1000nm/s). The start and end of each event in a curve were obtained through manual annotation. More statistical information on the data, including curve lengths and number of events per curve, is given in Table {#t:statistics}.

The expected rupture forces and loading rates were calculated from the annotated data. This process is shown graphically in Figure {#f:Rupture}. The region leading up to the event start was fit by a line, with the region length set to $\tau$ (see Table {#t:Parameters}). The loading rate was assigned to the line's slope, and the rupture force was calculated by the value of the line where the data in the region was last above the line. This is the same procedure used to calculate the loading rate and rupture force of predicted events in Figure {#f:Performance}. 

### Algorithm design and analysis

All timing and tuning results were obtained using a desktop with 16 GB of RAM, a 3.7 GHz i7 CPU, and a 1 TB hard drive. 

### Algorithm description

FEATHER improves on previous methods by using information present in the approach of the AFM cantilever to the surface-bound molecules (Figure {#f:FeatherExample}).  Figure {#f:Code} lists the pseudocode for the method. The algorithm is based on a probabilistic model of a signal lacking any events, called the *no-event model*, described in {#s:DesignDetails}. The algorithm has the following basic steps:

1. Estimate the no-event parameters (see Figure {#f:FeatherExample}) from the approach curve.
2. Fit the no-event model to the retract curve.
Calculate the upper bound on the probability of each retract point given the model.
3. Iteratively update the probability to remove false positives.
4. Report contiguous regions with probabilities lower than a user-specific threshold as events.

 FEATHER calculates an upper bound on the probability of no event occurring at each time point, using parameters estimated from the approach portion of the curve (see Section {#s:DesignDetails}) and a smoothing parameter from the user (see Table {#t:Parameters}). The probability at each point is iteratively updated to remove the effect of adhesions and other false positives. This last step is the only step requiring knowledge specific to SMFS. As shown in Figure {#f:FeatherExample}, the result is a probability at each time point which drops from one towards zero near events. A threshold probability is set by the user or optimized by a tuning routine (see Table {#t:Parameters} and Section {#s:Tuning}). Contiguous regions of time with probabilities below the threshold are considered having a single event, and the rupture properties are determined within each region as described in Section {#s:Annotation}.

Name       Meaning                            Value used in this work|
--------   --------------                     ----------------
$\tau$     Number of points for spline grid   2% of curve length
threshold  Probability used to reject events  Determined by Tuning (Figure {#f:Tuning})

[Table Caption]

### Choice of methods for comparison

The following two algorithms were chosen for comparison to FEATHER: 

- The AFM-specific 'event_find' routine from the OpenFovea AFM analysis package.@roduit_openfovea:_2012
- The general-purpose 'find_peaks_cwt' method from the Scientific Python package. @Jones_SciPy:_2001

 These methods were chosen to provide a representative sample of the viable techniques used in AFM data analysis. Unlike quadratic alignment algorithms in contour length space, these methods scale like O(N) and O(N$\cdot\log$(N)) respectively, where N is the length of a curve to be analyzed. Linear or near-linear scaling is desirable (see Figure {#f:Timing}) for analyzing hundreds of force-extension curves with millions of points each. These baselines represent two common approaches to event detection in AFM, since the OpenFovea method is a thresholding algorithm, and the Scientific Python method uses wavelet transformations. 

### Performance metrics

Two metrics were used for comparing the event-finding performance of FEATHER with the human-annotated data. The metrics reported are listed in Table {#t:metrics}. The event error metric, $P_{95}$, is the 95th percentile of relative error between predicted and true event locations (see Figure {#f:DistanceMetric}). The rupture Bhattacharya coefficient's complement, $1-<d_{(\nu,F),t}^{\frac{1}{2}}|d_{(\nu,F),p}^{\frac{1}{2}}>$, reports the mismatch between the true and predicted distribution over loading rates and rupture forces. 

### Algorithm tuning

All three algorithms were tuned using 5-fold cross validation. Cross validation was performed at fifteen log-spaced parameters over the useful parameter range of the algorithms (Figure {#f:Tuning}). The parameter value minimizing the Bhattacharya coefficient's complement for an algorithm was considered the algorithm's best parameter. Data shown in the results consists of the concatenation of all validation folds for each algorithm's best parameter. 


Since tuning the baselines on the full dataset would have required more than eight cpu-months (compared to $\approx$1.5 cpu-days for FEATHER, see Figure {#f:Timing}), a smaller subset of data was used for comparing the algorithms. In particular, the subset of the data with the smallest number of points per curve - 200 curves with v=1000}{nm/s}, N $\approx{}10^{5}$ (see Table {#t:statistics}) - was used for results comparing FEATHER to the baselines. FEATHER was also tuned separately on the larger, more complex dataset, with similar results to those reported in the rest of the paper (Figure {#f:LargeDataset}). This demonstrates that FEATHER generalizes well to a wide range of data sets sizes and experimental parameters.


![**Figure {#f:flowchart}.** Algorithm flowchart So that
blah blah blah
linebreak
blah blah](figures/flowchart.png)


## Results and Discussion

Table {#t:AppliedMetrics} lists the performance metrics for each algorithm. This section describes the metrics in more detail for each algorithm. 

Figure {#f:PerformanceFEATHER} shows the event detection performance of FEATHER. The distributions of the two types of relative distance errors, as defined in Figure {#f:DistanceMetric} and Table {#t:metrics}, exhibit a higher degree of symmetry than the baselines (see below). This symmetry is due to a low false positive rate, since false positives skew the distances from predicted to true points ($d_{\mathrm{p}\rightarrow\mathrm{t}}$). In addition, the error distribution is closer to zero than the baselines, indicating most errors are a small fraction of the length of the curve. Finally, Figure {#f:DistanceMetric} shows good agreement between the true two-dimensional distributions over rupture forces and loading rates and FEATHER's prediction.

The baselines do not perform as well as FEATHER. Figure {#f:PerformanceFovea} and Figure {#f:PerformanceScipy} are identical in structure to Figure {#f:PerformanceFEATHER}, but they describe the performance of the OpenFovea and Scientific Python baselines, respectively. OpenFovea's relative distance error distribution is heavily skewed towards one and asymmetric between the two error types, indicating many false positives. The same distribution for Scientific Python is more symmetric between error types, indicating fewer false positives than OpenFovea. However, the distance error distribution for Scientific Python has many points near 1, indicating a high error in predicted event locations. Neither baseline accurately reproduces the expected and predicted distribution over rupture forces and loading rates.

As defined by Table {#t:metrics}, relative to the best baseline FEATHER improves the relative and absolute event error by a factor of about XXX and improves the Bhattacharya coefficient's complement by a factor of about XXX. For completeness, Figure {#f:Performance} lists the performance of all three algorithms on a single page.


![**Figure {#f:flowchart}.** Algorithm flowchart So that
blah blah blah
linebreak
blah blah](figures/flowchart.png)

## Conclusion

 FEATHER provides order-of-magnitude improvements in physically-relevant metrics for algorithm performance relative to the baselines presented. As shown in Figure {#f:LargeDataset}, FEATHER also generalizes well to a wide range of SMFS experimental conditions.  

This work applies FEATHER to unfolding in SMFS data, but the algorithm could be generalized to a wider array of applications. A natural generalization to FEATHER would search for refolding events. Refolding experiments relax an unfolded molecule back to its native state by reducing the applied force. In addition, domain-specific models could be combined with FEATHER to further filter data, based on the specific experiment. Consider an experimenter attempting to fit a polymer model to the final rupture in the force-extension curve of a DNA hairpin of contour length 100}{nm}. After first obtaining all rupture locations in the data via FEATHER, the experimenter could fit a polymer model to each curve up to final rupture, obtain the contour length, and remove all data with contour lengths outside of (100$\pm \epsilon$)}{nm}, where $\epsilon$ is a small number. By predicting where events occur without providing any domain-specific model of the event, FEATHER is one tool in a longer SMFS analysis pipeline.

FEATHER could be used in different scientific domains than SMFS. To reduce the risk of overfitting to a specific type of SMFS experiment, FEATHER purposefully makes no assumptions about domain-specific data models (*e.g.* a polymer model for the force-extension relationship of DNA), except removing SMFS domain-specific false positives. The step which removes false positives could be changed or removed to fit the data application. SMFS data is a natural application for FEATHER, since each retract curve, with its scientifically interesting data, is paired with an approach curve which is assumed lacking any events. The parameters for FEATHER are estimated from the approach curve. FEATHER only requires estimating the relevant parameters (Figure {#f:FeatherExample}) and defining events by discontinuous time derivatives in otherwise continuous time series data. In scientific domains meeting these conditions, FEATHER could be a powerful tool for data analysis.  FEATHER's order-of-magnitude improvements in event detection improves the quality and reproducibility of event detection in piecewise-continuous time series data.

## Acknowledgments

## References

## Supplemental Material


![**Figure {#f:flowchart}.** Algorithm flowchart So that
blah blah blah
linebreak
blah blah](figures/flowchart.png)

A reference to Figure {#f:flowchart}. Here is a citation @dudko_theory_2008

