**Detecting Molecular Unbinding and Unfolding Events in Force Spectroscopy Data via a Bayesian Algorithm**

P. Heenan, R. Frongillo, J. Boyd-Graber, T. Perkins


ABSTRACT In single-molecule force spectroscopy (SMFS) experiments, mechanical forces are applied to individual biomolecules via a probe such as a bead or cantilever. Experiments which mechanically dissociate the secondary or tertiary structures of molecules are known as unfolding experiments. In unfolding experiments, a molecule binds to the probe, is pulled and possibly unfolded, and finally unbinds from the probe. The force on the molecule and its time derivative just before an unfolding or unbinding event are known as the event's rupture force and loading rate, respectively. The rupture force and loading rate at each event need to be known for common SMFS analyses. Identifying events in SMFS data is hindered by the presence of noise. This paper introduces a new algorithm, FEATHER (**F**orce **E**xtension **A**nalysis using a **T**estable **H**yptothesis for **E**vent **R**ecognition), to identify the locations of events in SMFS data. FEATHER features a relative event location error of XXX, a XXX-fold improvement relative to the best baseline, and a XXX for the distribution of loading rates and rupture forces of XXX, a factor of XXX improvement relative to the best baseline. As a linear-time algorithm for reproducible event identification, FEATHER improves the quality of analysis of SMFS data.  

Received for publication "Staff will complete"  and in final form "Staff will complete".

# Introduction {#sec:Intro}

In single-molecule force spectroscopy (SMFS) experiments, a force probe attaches to a molecule and stretches it while measuring force and extension over time (XXX ). These data are transformed into information such as kinetic rates of processive enzymes @comstock_direct_2015, protein-ligand bond strength @yuan_energy_2000, and the energy landscapes of proteins and nucleic acids @dudko_theory_2008. The location and properties of the *ruptures* in the data (see Figure {#fig:diagram}) are required for many common  analyses, such as applying polymer models and determining the molecular energy landscape.

----
![](figures/diagram.png){#fig:diagram}
**Figure** {#fig:diagram} Idiot Diagram
----


The analysis and interpretation of SMFS experiments is dependent on the attachment of the probe to the molecule of interest. Quality attachments between the probe and the sample may occur less than one curve in tens of thousands @bosshart_reference-free_2012. Strategies exist to improve attachment rates by coating the force probe with a molecule that binds readily to the sample of interest @walder_robert_rapid_nodate. However, the majority of data is still uninterpretable and must be removed. Excluding data without events and identifying the locations of events is a major challenge in SMFS data management and analysis.

Methods exist to automate detecting events within SMFS data. Automation removes the burden of manual data filtering and improves scientific reproducibility. Techniques for automated event detection in \fec{}s include aligning by the contour length at each extension by applying polymer models [@bosshart_reference-free_2012; @kuhn_automated_2005] thresholding based on signal or noise characteristics;[@gergely_semi-automatized_2001; @roduit_openfovea:_2012] and classification based on transformations of the data into frequency or derivative spaces [@kasas_fuzzy_2000; @garcia-masso_automated_2016; @benitez_searching_2017]. These methods do provide an increased degree of autonomy, but their use is limited by their lack of generalization. Thresholding or transmformation algorithms typicaly recquire optimizing many hard-to-interpret parameters, and contour-length alignment algorithms bias results towards dominant features and necessarily require a polymer model for the contour length (as a function of force and extension). 

This work describes a new method for detecting events in \fec{}s.  The algorithm, named FEATHER (**F**orce **E**xtension **A**nalysis using a **T**estable **H**yptothesis for **E**vent **R**ecognition), requires no \textit{a priori} knowledge of the polymer under study, does not bias data interpretation towards the dominant behavior of the data, and has two easy-to-interpret parameters which generalize well to typical SMFS data.

# Materials and Methods

--------------------------
![{#fig:flowchart}](figures/flowchart.png)
**Figure {#fig:flowchart}.** Algorithm flowchart
--------------------------

## Algorithm description

 FEATHER calculates an upper bound on the probability of no event occurring at each time point, using parameters estimated from the approach portion of the curve (see Section {#sec:DesignDetails}) and a smoothing parameter from the user (see Table {#tbl:Parameters}). The probability at each point is iteratively updated to remove the effect of adhesions and other false positives. This last step is the only step requiring knowledge specific to SMFS. As shown in Figure {#fig:flowchart}, the result is a probability at each time point which drops from one towards zero near events. A threshold probability is set by the user or optimized by a tuning routine (see Table {#tbl:Parameters} and Section {#sec:Tuning}). Contiguous regions of time with probabilities below the threshold are considered having a single event, and the rupture properties are determined within each region as described in Section {#sec:Annotation}.

## Choice of methods for comparison

The following two algorithms were chosen for comparison to FEATHER: 

- The thresholding 'event_find' routine from the OpenFovea AFM analysis package.@roduit_openfovea:_2012
- The wavelet-based 'find_peaks_cwt' method from Scientific Python. 
@jones_scipy:_2001

 These methods were chosen to provide a representative sample of the viable techniques used in AFM data analysis, since they utilize thresholding and wavelet transformation, two classes of event-detection techniques. 

## Performance metrics

Two metrics were used for comparing the event-finding performance of FEATHER with the human-annotated data. The metrics reported are listed in Table {#tbl:metrics}. The event error metric, $P_{95}$, is the 95th percentile of relative error between predicted and true event locations (see Figure {#fig:diagram}). The rupture Bhattacharya coefficient's complement reports the mismatch between the true and predicted distribution over loading rates and rupture forces.  (XXX label )


--------------------------------------------------- 
![{#fig:performance}](figures/performance.png)      
**Figure {#fig:performance}.** Performance figure   
--------------------------------------------------- 

# Results and Discussion

Table {#tbl:AppliedMetrics} lists the performance metrics for each algorithm. Figure {#fig:performance} shows the event detection performance of each algorithm. 

For FEATHER, the distributions of the two types of relative distance errors, as defined in Figure {#fig:diagram} and Table {#tbl:metrics}, exhibit a higher degree of symmetry than the baselines (see below). This symmetry is due to a low false positive rate, since false positives skew the distances from predicted to true points ($d_{\mathrm{p}\rightarrow\mathrm{t}}$). In addition, the error distribution is closer to zero than the baselines, indicating most errors are a small fraction of the length of the curve. Finally, Figure {#fig:performance} shows good agreement between the true two-dimensional distributions over rupture forces and loading rates and FEATHER's prediction.


The baselines do not perform as well as FEATHER. In Figure {#fig:performance}, OpenFovea's relative distance error distribution is heavily skewed towards one and asymmetric between the two error types, indicating many false positives. The same distribution for Scientific Python is more symmetric between error types, indicating fewer false positives than OpenFovea. However, the distance error distribution for Scientific Python has many points near 1, indicating a high error in predicted event locations. Neither baseline accurately reproduces the expected and predicted distribution over rupture forces and loading rates.

As defined by Table {#tbl:metrics}, relative to the best baseline FEATHER improves the relative and absolute event error by a factor of about XXX and improves the Bhattacharya coefficient's complement by a factor of about XXX. For completeness, Figure {#fig:performance} lists the performance of all three algorithms on a single page.

Name 	      	    | Rupture BCC ($\downarrow$) | Relative event error $P_{95}$ ($\downarrow$)
------------------- | ------------| --------------------------------
FEATHER             | **0.00501** | **0.00648**
OpenFovea 	    | 0.287 	  | 0.421
Scientific Python   | 0.0257 	  | 0.201
[{#tbl:Performance}]


# Conclusion

 FEATHER provides order-of-magnitude improvements in physically-relevant metrics for algorithm performance relative to the baselines presented. As shown in Figure XXX {#fig:LargeDataset}, FEATHER also generalizes well to a wide range of SMFS experimental conditions.  

This work applies FEATHER to unfolding in SMFS data, but the algorithm could be generalized to a wider array of applications. A natural generalization to FEATHER would search for refolding events. Refolding experiments relax an unfolded molecule back to its native state by reducing the applied force. In addition, domain-specific models could be combined with FEATHER to further filter data, based on the specific experiment. Consider an experimenter attempting to fit a polymer model to the final rupture in the force-extension curve of a DNA hairpin of contour length 100}{nm}. After first obtaining all rupture locations in the data via FEATHER, the experimenter could fit a polymer model to each curve up to final rupture, obtain the contour length, and remove all data with contour lengths outside of (100$\pm \epsilon$)}{nm}, where $\epsilon$ is a small number. By predicting where events occur without providing any domain-specific model of the event, FEATHER is one tool in a longer SMFS analysis pipeline.

FEATHER could be used in different scientific domains than SMFS. To reduce the risk of overfitting to a specific type of SMFS experiment, FEATHER purposefully makes no assumptions about domain-specific data models (*e.g.* a polymer model for the force-extension relationship of DNA), except removing SMFS domain-specific false positives. FEATHER only requires estimating the relevant parameters and defining events by discontinuous time derivatives in otherwise continuous time series data. In scientific domains meeting these conditions, FEATHER could be a powerful tool for data analysis.  FEATHER's order-of-magnitude improvements in event detection improves the quality and reproducibility of event detection in piecewise-continuous time series data.

# Acknowledgments

# References

[//]: # (see: 
stackoverflow.com/questions/16427637/pandoc-insert-appendix-after-bibliography)

<div id="refs"></div>

# Supplemental Information

\newpage