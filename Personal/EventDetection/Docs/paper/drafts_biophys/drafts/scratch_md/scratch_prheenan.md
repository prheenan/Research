**Detecting Molecular Unbinding and Unfolding Events in Force Spectroscopy Data via a Bayesian Algorithm**

P. Heenan, R. Frongillo, J. Boyd-Graber, T. Perkins


ABSTRACT In single-molecule force spectroscopy (SMFS) experiments, mechanical forces are applied to individual biomolecules via a probe such as a bead or cantilever. Experiments which mechanically dissociate the secondary or tertiary structures of molecules are known as unfolding experiments. In unfolding experiments, a molecule binds to the probe, is pulled and possibly unfolded, and finally unbinds from the probe. The force on the molecule and its time derivative just before an unfolding or unbinding event are known as the event's rupture force and loading rate, respectively. The rupture force and loading rate at each event need to be known for common SMFS analyses. Identifying events in SMFS data is hindered by the presence of noise. This paper introduces a new algorithm, FEATHER (**F**orce **E**xtension **A**nalysis using a **T**estable **H**yptothesis for **E**vent **R**ecognition), to identify the locations of events in SMFS data. FEATHER features a relative event location error of 0.006, a 30-fold improvement relative to the best baseline, and a Bhattacharya coefficient's complement for the distribution of loading rates and rupture forces of 0.005, a factor of 5 improvement relative to the best baseline. As a linear-time algorithm for reproducible event identification, FEATHER improves the quality of analysis of SMFS data.  


# Introduction {#label_sec:Intro}


In single-molecule force spectroscopy (SMFS) experiments, a force probe attaches to a molecule and stretches it while measuring force and extension over time, known as a force-extension curve (Figure {#ref_fig:diagram}). These data are transformed into information such as kinetic rates of processive enzymes @comstock_direct_2015, protein-ligand bond strength @yuan_energy_2000, and the energy landscapes of proteins and nucleic acids @dudko_theory_2008. The location and properties of the *ruptures* in the data (Figure {#ref_fig:diagram}) are required for many common  analyses, such as applying polymer models and determining the molecular energy landscape.


----
![](figures/diagram.png){#label_fig:diagram}
**Figure {#ref_fig:diagram}** Energy landscape characterization of protein unfolding.  **(A)** A cartoon of a single-barrier energy potential and the polyprotein construct used in this paper. **(B)** An illustrative force versus time curve for the cartoon shown in (A). **(C)** An cartoon distribution of rupture forces, $F_\mathbf{R}$, with a model fit overlayed @dudko_theory_2008. **(D)** Applying the model from (C) to rupture force and loading rate yields the energetic parameters in (A). **(E)** An experimental force versus time curve used in this study, with events marked by green arrows. **(F)** A detailed subplot of a single event from (E), defining $d_{t\rightarrow p}$, the distance from a 'true' (annotated) event to the closest predicted event, and $d_{p\rightarrow t}$, the distance from a predicted event to the closest true event. Predicted events are for illustrative purposes only.
----


Methods exist to automate detecting events within SMFS data. Automation removes the burden of manual data filtering and improves scientific reproducibility. Techniques for automated event detection in SMFS data include aligning by the contour length at each extension by applying polymer models [@bosshart_reference-free_2012; @kuhn_automated_2005] thresholding based on signal or noise characteristics;[@gergely_semi-automatized_2001; @roduit_openfovea:_2012] and classification based on transformations of the data into frequency or derivative spaces [@kasas_fuzzy_2000; @garcia-masso_automated_2016; @benitez_searching_2017]. These methods do provide an increased degree of autonomy, but their use is limited by their lack of generalization. Thresholding or transmformation algorithms typicaly recquire optimizing many hard-to-interpret parameters, and contour-length alignment algorithms bias results towards dominant features and necessarily require a polymer model for the contour length (as a function of force and extension). 

This work describes a new method for detecting events in force-extension curves.  The algorithm, named FEATHER (**F**orce **E**xtension **A**nalysis using a **T**estable **H**yptothesis for **E**vent **R**ecognition), requires no \textit{a priori} knowledge of the polymer under study, does not bias data interpretation towards the dominant behavior of the data, and has two easy-to-interpret parameters which generalize well to typical SMFS data.

# Materials and Methods

## Data used to determine performance {#label_sec:datasets}

The following two data sets were hand-annotated for the purposes of determining algorithm event detection performance:

- A publically available polyprotein dataset (XXX)
- A 650nm dsDNA dataset taken for this paper (see Section {#ref_sec:SampleDetails}). 

Statistics on the data sets and the annotated events are in Section {#ref_sec:SampleDetails}

## Description of FEATHER

 FEATHER calculates an upper bound on the probability of no event occurring at each time point, using a smoothing parameter $\tau\in[0,1]$ from the user (see Table {#ref_tbl:Parameters}), a threshold probability defining the minimum probability a non-event can have, and parameters determined from the force versus time curve (see Figure S{#ref_fig:algorithm_details}). This process is illustrated in Figure {#ref_fig:flowchart} and described in more detail in Section S{#sec:DesignDetails}.


--------------------------
![{#label_fig:flowchart}](figures/flowchart.png)
**Figure {#ref_fig:flowchart}.** FEATHER’s algorithmic identification of rupture events in force. **(A)** A force versus time curve with a spline fit overlayed. **(B)** The probability of no event obtained by applying Chebychev's inequality to (A), as described in the text. **(C)** Transforming (B) to remove regions near the surface or with positive force derivatives. **(D)** Transforming (C) to remove regions where the force change is negligible, as described in {#sec:DesignDetails}. **(E-H)** The events and magnified regions obtained from each region less than a user-specified threshold (D). Plotting conventions are as in Figure {#ref_fig:diagram}.
--------------------------


## Choice of methods for comparison

The following algorithms were chosen for comparison to FEATHER: 

- The thresholding 'event_find' routine from the OpenFovea AFM analysis package.@roduit_openfovea:_2012
- The wavelet-based 'find_peaks_cwt' method from Scientific Python. 
@jones_scipy:_2001

 These methods were chosen to provide a representative sample of the viable techniques used in AFM data analysis, since they respectively utilize thresholding and wavelet transformation, two classes of event-detection techniques. 

## Performance metrics

Two metrics were used for comparing the event-finding performance of FEATHER with the human-annotated data. The metrics reported are listed in Table {#tbl:metrics}. The event error metric, $P_{95}$, is the 95th percentile of relative error between predicted and true event locations (see Figure {#ref_fig:diagram}). The rupture Bhattacharya coefficient's complement reports the mismatch between the true and predicted distribution over loading rates and rupture forces.


--------------------------------------------------- 
![{#label_fig:performance}](figures/performance.png)      
**Figure {#ref_fig:performance}.** Figure reporting FEATHER's 10x performance gain. **(A)** The histogram of relative errors between expected and predicted event locations for FEATHER, Open Fovea, and Scientific Python. Distributions skewed towards one indicate high error. The dotted line gives the 95% percentile of error. **(B)** The histogram of predicted and expected loading rates for FEATHER, Open Fovea, and Scienctific Python. **(C)** As (B), but for rupture forces. All expected values are from human-annotated data.
--------------------------------------------------- 

# Results and Discussion

Table {#ref_tbl:AppliedMetrics} lists the performance metrics for each algorithm on he polyprotein dataset. Figure {#ref_fig:performance} shows the event detection performance of each algorithm. As defined by Table {#tbl:metrics}, relative to the best baseline FEATHER improves the relative and absolute event error by a factor of about 30 and improves the Bhattacharya coefficient's complement by a factor of about 5. Additional tests on the dsDNA dataset (See Section {#ref_sec:datasets} and Figure {#ref_fig:DNA}) show even greater performance gains. FEATHER's performance detecting events in polyprotein and dsDNA data demonstrate FEATHER's general applicability and order-of-magnitude improvements relative to the baselines presented. 


Name 	      	    | Rupture BCC ($\downarrow$) | Relative event error $P_{95}$ ($\downarrow$)
------------------- | ------------| --------------------------------
FEATHER             | **0.00501** | **0.00648**
OpenFovea 	    | 0.287 	  | 0.421
Scientific Python   | 0.0257 	  | 0.201
[{#label_tbl:Performance}]


# Conclusion

This work applies FEATHER to unfolding in SMFS data, but the algorithm could be generalized to a wider array of applications. A natural generalization to FEATHER would search for refolding events. Refolding experiments relax an unfolded molecule back to its native state by reducing the applied force. In addition, domain-specific models could be combined with FEATHER to further filter data, based on the specific experiment. For example, after first obtaining all rupture locations in a dataset via FEATHER, an experimenter could fit a polymer model to each rupture and filter the results based on expected contour lengths. By predicting where events occur without providing any domain-specific model of the event, FEATHER is one tool in a longer SMFS analysis pipeline.

FEATHER could be used in different scientific domains than SMFS. To reduce the risk of overfitting to a specific type of SMFS experiment, FEATHER purposefully makes no assumptions about domain-specific data models (*e.g.* a polymer model for the force-extension relationship of DNA), except removing SMFS domain-specific false positives. FEATHER only requires estimating the relevant parameters and defining events by discontinuous time derivatives in otherwise continuous time series data. In scientific domains meeting these conditions, FEATHER could be a powerful tool for data analysis.  FEATHER's order-of-magnitude improvements in event detection improves the quality and reproducibility of event detection in piecewise-continuous time series data.

# Location and accessibility

FEATHER is written in python, with interfaces written for Matlab and Igor Pro. The code is available with working examples at <DOI>.

# References

[//]: # (see: 
stackoverflow.com/questions/16427637/pandoc-insert-appendix-after-bibliography)

<div id="refs"></div>

# Supplemental Information

\newpage