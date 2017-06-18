

# {#label_sec:SampleDetails} Cantilever and dsDNA sample preparation 

Site-specific chemistry is used to improve the acquisition rate and quality of data. The procedure for surface, sample, and cantilever preparation is described briefly in Appendix A and in detail elsewhere @walder_robert_rapid_nodate. Briefly, through polymerase chain reaction, double-stranded DNA is functionalized with dibenzocyclooctyl (DBCO) at the 5' end of one DNA strand to ensure a covalent bond with an azide-functionalized surface. The DNA is also functionalized with biotin at the other 5' end to ensure a specific but reversible bond with a streptavidin-coated cantilever. These two bonds ensure the tip-DNA bond is broken before the surface-DNA bond, preventing tip contamination. 

## Atomic force microscopy

All atomic force microscopy (AFM) measurements were carried out using an Asylum Cypher (Asylum Research Cypher ES). The spring constant and sensitivity were determined using the equipartition theorem method @hutter_calibration_1993. All measurements were carried out with a 2 second pause at the surface. To promote molecular attachments to the tip, measurements were taken over a series of square, 30 $\mathrm{\mu}$m grids at points separated by 1 $\mathrm{\mu}$m. The stage was translated to a new grid point after measuring three force-extension curves. 


## Azide-functionalized surfaces

Glass surfaces were etched using potassium hydroxide. 12 mm diameter glass disks (Ted Pella 26023) were sonicated (Branson 2510) for 5 minutes in 250 mL of ACS grade acetone (Fisher A18-4), 5 minutes in 250 mL of 95\% ethanol (Decon 2801), and transferred to a solution of 80 g of potassium hydroxide (Fisher P250-500) dissolved in  170 mL of 95\% ethanol and 80 mL of deionized water (Barnstead GenPure Pro) for 3 minutes. To removed residual KOH, the glass was serially diluted in two 1 L beakers of 18.2 M$\Omega$ deionized water. The etched surfaces were dried with 99.8\% pure nitrogen gas (Airgas NI UHP-300) and stored at room temperature in a dust-proof box. 

To create azide-functionalized surfaces, etched glass was activated using 30 minutes exposure in an UV-ozone chamber (Novascan PSDP Digital UV Ozone System) and incubated for three hours in a 70 mL, $60^{\circ}$C solution of 0.15 mg/mL 600 Dalton (Silane-PEG-Azide Nanocs PG2-AZSL-600) dissolved in toluene (Sigma 179418-4L). The surfaces were then mounted in a custom-made teflon-holder and rinsed serially in 250 mL of toluene, isopropanol (Fisher A416-4), and deionized water. The azide-functionalized surfaces were dried with nitrogen gas and stored in a dust-proof container at $4^{\circ}$C.  


Name           | Sequence
-------------- | -------------------------------- 
Forward primer | Dagttgttcctttctattctcactccgc     
Reverse primer | BtcaataaTcggctgtctTtccttatcaTtc  
[{#label_tbl:Sequences}Sequences for single-stranded DNA construct and for double-stranded DNA primers. The forward and reverse primers amplify a 647 nm piece of DNA, from positions 1607 to 3520 on the M13mp18 plasmid, as discussed in the text. Unmodified DNA bases are lowercase. The uppercase letters 'T','B', and 'D' respectively represent biotinylated Thymidine residues, a terminal Biotin separated from the sequence by triethyleneglycol spacer, and a terminal dibenzocyclooctyl (DBCO) separated from the sequence by triethyleneglycol spacer.]

## {#label_sec:Sample} DNA Samples

For the 1914 basepair (bp) double-stranded DNA, the 1607 forward-sense and 3520 reverse-sense primers (Table S{#ref_tbl:Sequences}) for the M13mp18 plasmid (New England BioLabs N4018S) were obtained from Integrated DNA Technologies  and used for polymerase-chain reaction (Millipore 71086-4) using 40 cycles on a thermocycler (Bio-Rad T100). The reverse-sense primer was modified to include three biotinylated thymidine bases and a 5-prime biotin after a PEG spacer. The forward-sense primer was modified to include a 5-prime dibenzocyclooctyl (DBCO) after a PEG-spacer (see Table S{#ref_tbl:Sequences}). After the PCR product was purified (Qiagen 28106) and the 1.9 kbp band selected by 2\% gel electrophoresis (Sigma A9539-500) as shown in Figure S{#ref_fig:Prep}, the agarose was eluted out (Bio-Rad 732-6165), the DNA solution was concentrated (Millipore UFC501096) and purified (Qiagen 28106). The purified construct was stored at $4^{\circ}$C in a solution of 10 mM Tris (Fisher BP151-1) and 1 mM EDTA (Fisher S311-500) at pH 8.0, referred to after as 'TE'. The typical recovery efficiency of DNA purification after polymerase chain reaction was 25\%. 

The purity of the DNA was verified by depositing 20 pmol of DNA in imaging buffer \textemdash{} 3 mM Nickel (Sigma 654597), 10 mM Hepes (Sigma H4034) at pH 7 \textemdash{} onto freshly-cleaved, 10 mm diameter, V1 mica (Ted Pella 50) for 10 minutes, rinsing serially 5 times with 1 mL of deionized water and 5 times with 1 mL of imaging buffer. The sample was then imaged with the Cypher AFM using a cantilever with a nominal 2 nm radius, a spring constant of $\approx 300 \frac{\text{pN}}{\text{nm}}$ Bruker SNL-10, a line scan rate of 1 Hz, a 1 \textmu{m} scan size, and free amplitude of $\approx 1$ nm (Figure S{#ref_fig:Prep}). 

For DNA deposition onto an azide surface, 20 \textmu{L} of the DNA at 40 nM was mixed with 80 uL of TE and deposited onto azide-functionalized glass affixed with 5 minute epoxy (Devcon 14250) to specimen disks (Ted Pella 16218) and incubated at $4^{\circ}$C overnight. The surfaces were rinsed by 7 mL of TE at pH 8.0 and 7 mL of Phosphate Buffered Saline (Millipore 524650), referred to after as PBS, with 1 mM EDTA at pH 7.4, and stored at $4^{\circ}$C. 

## Streptravidin-functionalized cantilevers

Functionalized cantilevers with nominal spring constants of 4 $\frac{\text{pN}}{nm}$  were used for all experiments. The protocol below etches away gold to improve force stability @sullan_atomic_2013 and covalently attaches streptatividin to the etched cantilevers to improve attachment to the biotinylated DNA. AFM cantilevers (Asylum Research BL-RC-150VB) were serially rinsed for 30 seconds in 50 mL of gold etchant (Transene TFA), 250 mL of deionized water, 50 mL of chromium etchant (Transene 1020), and 250 mL of deionized water. The tips were then serially rinsed for 30 seconds in 50 mL of deionized water, isopropanol, toluene, isopropanol, and deionized water again. After drying, the tips were activated via 30 minutes exposure in a UV-Ozone chamber and incubated for three hours in a 70 mL, $60^{\circ}$C solution of 0.15 mg/mL 600 Dalton Silane-PEG-Maleimide (Nanocs PG2-MLSL-600) dissolved in toluene. The maleimide-functionalized tips were serially rinsed in 50 mL of toluene, isopropanol, and water, immediately dried (KimTech 5511), and immersed for three hours in a 0.2 mg/mL solution of thiol-streptavidin (Protein Mods SAVT) in PBS at pH 6.75 with 1 mM TCEP (Thermo Scientific 77720) at room temperature in a moisture-proof container. Subsequently, the tips were transferred to $4^{\circ}$C overnight. After the $4^{\circ}$C incubation, to remove free streptavidin, the tips were serially rinsed in two 10 mL beakers of PBS at pH 7.4 and submerged in a 20 mL petri dish of PBS for 10 minutes. These functionalized tips were stored in 50 \textmu{L} of PBS at $4^{\circ}$C in plastic wafers (Entegris H22-100/101-0615) until loading into the atomic force microscope. 



# {#label_sec:DesignDetails}  Algorithm design

This section details the event-detection algorithm. The following conventions are followed:

\begin{enumerate}
 \item All variables with uppercase letters are random variables.
 \item All variables with lowercase letters are either instances of the corresponding uppercase random variables (\textit{i.e.}, measurements) or pure functions.
 \item All variables with a hat (\textit{e.g.}, $\hat{\epsilon}_t$) are estimators of a variable.
\end{enumerate}

--------------------------------------------------- 
![](figures/algorithm.pdf.png)      
{#label_fig:algorithm_details} Demonstrating how FEATHER works. **(A,B)** The approach and retract forces versus time, with spline fits overlaid. **(C,D)** The approach and retract $r_t$ versus time, with $\Sigma_t$ overlaid, demonstrating the signal-to-noise benefit of using $\Sigma_t$. **(E,F)** The approach and retract $\Sigma_t$, with the estimates of $\epsilon$ and $\sigma$ from the approach overlaid. **(G)** The probability of no event at each time has a sharp minima near the expected event location. For all subplots, the retract changes color at the location of a tagged event. 
--------------------------------------------------- 

FEATHER improves on previous methods by using information present in the approach of the probe to the surface-bound molecules (Figure S{#ref_fig:algorithm_details}). The algorithm is based on a probabilistic model of a signal lacking any events, called the *no-event model*, described in Section S{#ref_sec:DesignDetails}. The algorithm has the following basic steps:


1. Estimate the no-event parameters (see Figure S{#ref_fig:algorithm_details}) from the approach curve.
2. Fit the no-event model to the retract curve.
3. Calculate the upper bound on the probability of each retract point given the model.
4. Iteratively update the probability to remove false positives.
5. Report contiguous regions with probabilities lower than a user-specific threshold as events.

FEATHER uses a probabilistic model for the portion of the force-extension curve where force is applied to the molecule of interest, referred to as the 'retract', based on the portion of the force-extension curve when the probe is not in contact with the molecule, referred to as the 'approach'. The algorithm fits and subtracts a smoothing spline from the approach, yielding an expected mean and variance of the residual's standard deviation within a window of $\pm\tau$. Applying this procedure to the retract yields a residual mean standard deviation at each point in time. This residual is transformed into a probability using Chebyshev's inequality and the expected mean and variance from the approach (Equation S{#ref_eq:featherprobability}). This probability at each point is iteratively updated to remove the effect of adhesions and other false positives. As shown in Figure {#ref_fig:flowchart}, the result is a probability at each time point which drops from one towards zero near events. A threshold probability is set by the user or optimized by a tuning routine (see Table S{#ref_tbl:Parameters} and Section S{#ref_sec:Tuning}). 

Contiguous regions of time with probabilities below the threshold are considered a single event, and the rupture properties are determined within each region. In the region, let the first time the smoothing spline's derivative crosses its median within the region is defined as index 'i'. Then, the region leading up to index i was fit by a line, with the region length set to $\tau$ (see S{#ref_tbl:Parameters}).  The loading rate was assigned to the line's slope, and the rupture force was calculated by the value of the line where the data in the region was last above the loading rate line.

##Defining the no-event hypothesis

FEATHER defines an event as a discontinuity in piecewise-continuous time series data. In SMFS, events occur when the force applied to a molecule exhibits a discontinuity as the molecule passes over an energy barrier. FEATHER assumes that events take place on time scales much smaller than the response time of the probe. If this is not true, the data is assumed smoothed until this condition is reached. The algorithm models the data assuming no event is occurring at a given time '$t$' and finds where the probability of a measurement is low. Hereafter, the definition of an event and these assumptions will be referred to as the *no-event hypothesis*. 

##Mathematical background for testing the no-event hypothesis

Under the no-event hypothesis, the noise-dependent distribution of force '$F_t$' for a discrete series of forces sampled at time points 't' can be well-approximated by the sum of a smooth signal and a noise distribution (Figure S{#ref_fig:algorithm_details}):

$F_t = g_t + X(0,\sigma^2)$

where $g_t$ is a smooth function with a continuous first derivative and X is a random variable with zero mean and variance $\sigma^2$. The closed form of the smooth signal and the noise distribution are assumed unknown and can vary from one \fec{} to the next. If '$g^{*}_t$' is a function with a continuous first derivative approximately equal to $g_t$ such that $\forall t,\epsilon_t\equiv|g^{*}_t-g_t|$ for real $\epsilon_t\ge 0$, then the error distribution $R_t$ is defined such that: 

$R_t \equiv F_t - g^{*}_t = \epsilon_t + X(0,\sigma)$

where

$E[R_t^2] -E[R_t]^2 = [(g_t-g^{*}_t)^2 + E[X(0,\sigma)^2]] - (g_t-g^{*}_t)^2  = \sigma^2$ {#label_eq:feathersigma}

and 

$|E[R_t]| = \epsilon_t \le E[|R_t|]$   {#label_eq:featherepsilon}

Under these assumptions, the probability 'P' of measuring $r_t$ is bounded by Chebyshev's inequality:

$P( |R_t-\epsilon_t| \ge |r_t-\epsilon_t| ) \le
(\frac{\sigma}{|r_t-\epsilon_t|})^2$ {#label_eq:featherprobability}


For \emph{any} noise distribution, Equation S{#ref_eq:featherprobability} bounds the probability of a measurement under the no-event hypothesis, given the force approximation $g^{*}_t$ (which in turn yields the estimator error $\hat{\epsilon}$ and the noise $\hat{\sigma}$ by Equation S{#ref_eq:featherepsilon} and Equation {#eq:feathersigma}). A low probability at a given time implies the measurement is unlikely under the no-event hypothesis. 

##Accurate estimators for hypothesis testing

To obtain an accurate estimator for $g_t$, the data must be smoothed. The approximation to the noiseless force $g^{*}_t$ is obtained by fitting a least-squares second-order basis spline @dierckx_algorithm_1975 to the force versus time curve. The spline is second-order to ensure a continuous first derivative (Figure S{#ref_fig:algorithm_details}), and the spline knots are spaced uniformly at intervals of  the user-specified $\tau$ (see Table S{#ref_tbl:Parameters}). Figure S{#ref_fig:algorithm_details} is a representative demonstration of the spline fitting. Determining  $g^{*}_t$ immediately gives $\hat{r}_t$.

Using $r_t$ as shown in Equation S{#ref_eq:featherepsilon} does not provide a strong signal in the presence of an event (see Figure S{#ref_fig:algorithm_details}). In order to improve the method, $r_t$ was replaced by the distribution of windowed standard deviations '$\Sigma$'. $\Sigma$ is defined as the standard deviation of $r_t$ centered at t with a window of $[-\frac{\tau}{4},\frac{\tau}{4}]$. Using $\Sigma$ instead of $r_t$ provides a much stronger signal in the presence of an event (see Figure S{#ref_fig:algorithm_details}).  

The noise variables $\sigma$ and $\epsilon$ are estimated from the distribution of standard deviations $\Sigma$ on the region of the approach curve where the AFM tip is not in contact with the surface. From this distribution, $\hat{\sigma}$ is set to the standard deviation of $\Sigma$, and $\hat{\epsilon}_t$ is approximated by the median. The median is used instead of the mean to remove the influence of possible false positive events in the approach. The removal of these pseudo-events is necessary to ensure accurate estimators for $\sigma$ and $\epsilon$, which are based on the no-event hypothesis. 

The quality of FEATHER's results are improved by multiplying the no-event probability, as discussed above, by the integral force, force derivative, and force differential Chebyshev probabilities. The calculation of each of these probabilities is exactly the same as Equation {#eq:featherprobability}, with the variables changed appropriately. Specifically, the relevant operation (integration, differentiation, or force difference) is applied to the approach, estimates for the operation-specific $\epsilon$ and $\sigma$ are obtained, yielding the operation-specific probability.



--------------------------------------------------- 
![](figures/prep.pdf.png)      
{#label_fig:prep} Purity of the sample preparation. **(A)** 2\% electrophoretic agarose gel, showing a major band just below 2kbp, as expected. **(B)** A false-color AFM image of the DNA bound to mica. 
--------------------------------------------------- 


## Data Annotation

Two hundred force-extension curves with events were obtained at three pulling velocities (100 nm/s, 500 nm/s, 1000nm/s). The start and end of each event in a curve were obtained through manual annotation. More statistical information on the data, including curve lengths and number of events per curve, is given in Table S{#ref_tbl:statistics}.


v (nm/s) | $N_\mathrm{curves}$ | $\mu_{\mathrm{Curve Size}}$ | $\sigma_{\mathrm{Curve Size}}$ | $N_{\mathrm{e}= 1}$ | $N_{\mathrm{e}= 2}$ | $N_{\mathrm{e}= 3}$ | $N_{\mathrm{e}\ge4}$ 
----- | --- | ------ | ----- | --- | -- | - | - |
100 | 200 | 666000 | 46000 | 160 | 31 | 5 | 4
500 | 200 | 200000 | 1000 | 142 | 39 | 10 | 9
1000 | 200 | 116000 | 7000 | 176 | 22 | 2 | 0
[{#label_tbl:statistics} Statistical information on the 650nm dsDNA data set. For each loading rate v in the data set, this table lists the number of curves $N_{\mathrm{curves}}$; mean and standard deviation of curve sizes, $\mu_{\mathrm{Curve Size}}$ and $\sigma_{\mathrm{Curve Size}}$, in data points; the number of curves with `x' events $N_{\mathrm{e=x}}$ for x$\in\{1,2,3\}$; and the number of curves with greater than or equal to 4 events, $N_{\mathrm{e}\ge4}$. ]



v [nm/s] | $N_\mathrm{curves}$ | $\mu_{\mathrm{Curve Size}}$ | $\sigma_{\mathrm{Curve Size}}$ | $N_{\mathrm{e}= 3}$ | $N_{\mathrm{e}= 4}$ | $N_{\mathrm{e}= 5}$ | $N_{\mathrm{e}= 6}$ | $N_{\mathrm{e}\ge7}$
----- | --- | ------ | ----- | --- | -- | - | - | - |
Various | 152 | 93000 | 26000 | 7 | 24 | 48 | 66 | 7
[{#label_tbl:protein_statistics} Statistical information on the polyprotein data set. Conventions are as in Table S{#ref_tbl:statistics}.]



--------------------------------------------------- 
![](figures/landscape.pdf.png)
{#label_fig:DNA} On the dsDNA dataset, FEATHER has orders-of-magnitude better performance compared to the baseline algorithms. **(A1)** The distribution of distances from predicted to true points, $d_{\mathrm{p}\rightarrow\mathrm{t}}$, and from true to predicted points, $d_{\mathrm{t}\rightarrow\mathrm{p}}$, for FEATHER. **(A2)** FEATHER's two-dimensional distribution of true and predicted rupture forces and loading rates, as defined in Table S(#ref_tbl:metrics). The range of the plot is limited to the middle 98 percent of the data. **(A3,A4)** The histograms of rupture forces and loading rates, respectively, for FEATHER. The range of these plots are limited as in (B). **(A5)** The metrics defined in Table S(#ref_tbl:metrics) applied to FEATHER.. **(B1-B5)** As A1-A5, except for the Open Fovea baseline. **(C1-C5)** As A1-A5, except for the Scientific Python baseline.
--------------------------------------------------- 

# {#label_sec:tuning} Algorithm tuning

All three algorithms were tuned using 5-fold cross validation. Cross validation was performed at fifteen log-spaced parameters over the useful parameter range of the algorithms. The parameter value minimizing the Bhattacharya coefficient's complement for an algorithm was considered the algorithm's best parameter. Data shown in the results consists of the concatenation of all validation folds for each algorithm's best parameter. 


--------------------------------------------------- 
![](figures/FEATHER_full.pdf.png)      
{#label_fig:full_dataset} FEATHER generalizes well to a wide range of Data. The subplots in this figure are formatted as in Figure S{#ref_fig:DNA} but with FEATHER applied to the full data set listed in Table S{#label_tbl:statistics}. 
--------------------------------------------------- 

Since tuning the baselines on the full dataset would have required more than eight cpu-months (compared to $\approx 1.5$ cpu-days for FEATHER, see Figure S{#ref_fig:timing}), a smaller subset of data was used for comparing the algorithms. In particular, the subset of the data with the smallest number of points per curve - 200 curves with v=1000}{nm/s}, N $\approx{}10^{5}$ (see Table S{#ref_tbl:statistics}) - was used for results comparing FEATHER to the baselines. FEATHER was also tuned separately on the larger, more complex dataset, with similar results to those reported in the rest of the paper (Figure S{#ref_fig:full_dataset}). This demonstrates that FEATHER generalizes well to a wide range of data sets sizes and experimental parameters.


Name       Meaning                            Value used in this work
--------   --------------                     ----------------
$\tau$     Number of points for spline grid   2% of curve length
threshold  Probability used to reject events  Determined by parameter sweep


# Performance metrics


Name             | Notation			           | Range | Optimum 
---------------- |--------------------------------------   | ---   | ---
K     | number of curves				| -     | -
$N_k$ | points in curve k				| -     | -
$d_{i\rightarrow j,k}$ | distribution of pointwise distances in 'k' from 'i' ruptures to the closest 'j' rupture or $N_k$ if none    	      	| -      | -
$P_{x}$ |  the 'x'-th percentile of the concatenation of $\frac{1}{N_k}d_{t\rightarrow p,k}$ and $\frac{1}{N_k}d_{p\rightarrow t,k}$ over all k | -      |  -
$\nu_i$ | histogram of 'i' loading rates over all k     | -	 | - 
$F_i$ | histogram of 'i' rupture forces over all k      | -	 | - 
$d_{(\nu,F),i}$ | joint histogram of $\nu_i$ and $F_i$ divided by K  | -      | -
**relative event error** | $P_{95}$     	        | [0,1]  | 0 
**rupture BCC** | 1-$<d_{(\nu,F),t}^{\frac{1}{2}}|d_{(\nu,F),p}^{\frac{1}{2}}>$ | [0,1] | 0 
[{#label_tbl:metrics} The definitions of the performance metrics reported. The metrics are bolded, and the quantities that they depend on are listed first. BCC stands for Bhattacharya coefficient's complement.  Throughout, 'k' refers to the index of a force-extension curve, and 'i' and 'j' refer to either true or predicted. For example, $d_{t\rightarrow p,4}$ represents the distances from the true to the predicted events in force-extension curve 4. ]


# {#label_sec:Timing} Algorithmic time complexity

All timing and tuning results were obtained using a desktop with 16 GB of RAM, a 3.7 GHz i7 CPU, and a 1 TB hard drive. 

Figure S{#ref_fig:timing} compares the runtimes, T(N), of FEATHER and the baselines. The runtime of each algorithm is linear with curve size. FEATHER has a roughly tenfold better asymptotic slope than the baselines.  

--------------------------------------------------- 
![](figures/timing.pdf.png)      
{#label_fig:timing} Big-O runtime of each algorithm is approximately linear. **(A)** The runtime per curve versus number of points in the curve, T(N), is reported for each algorithm. Lines around each algorithm show upper and lower bounds of the form $T_{\mathrm{upper}}$(n) = $a_0$ + 2$a_1$n  and $T_{\mathrm{lower}}$(n) = $a_0$ + $\frac{1}{2}a_1$n, where $a_0$ and $a_1$ are coefficients to the least-squares fit of the runtime. **(B)** This plot shows the asymptotic runtime per data point, $a_1$, for each algorithm from (A).
--------------------------------------------------- 

