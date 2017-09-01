**Free Energy Landscape of Bacteriorhodopsin Reveals Local Variations in Unfolding Energy**

P. Heenan, H. Yu, M. Siewny , T. Perkins

ABSTRCT 


Determining structural causes for membrane protein stability is a historically difficult problem. Single molecule force spectroscopy (SMFS) experiments can mechanically pull apart a molecule and determine its associated free energy landscape, and hence stability, as a function of molecular extension. Atomic force microscopy (AFM) is one SMFS technique that unfolds membrane proteins in their native lipid bilayers. Traditionally, AFM-SMFS data has been limited by poor instrumental time resolution, making state identification difficult, and the perturbing effect of the AFM cantilever, causing systematic overestimation of the landscape. Recent advancements in ultra-short ($<$L$>$=?m) cantilevers with 1$\upmu$s time resolution can detect previously obscured structural states within SMFS data. In addition, developments in data analysis techniques such as the inverse Weierstrass transform approximately remove the perturbation of the force probe on the landscape. By combining the advances in algorithmic analysis with improvements in cantilever technology, we reconstructed the energy landscape of the membrane protein bacteriorhodopsin as a function of molecular end-to-end distance with 2 angstrom resolution. The average energy of unfolding per amino acid was 0.7 $\pm$ 0.1 kcal/mol, in agreement with past work. However, there was significant variation in unfolding energy cost per amino acid within bacteriorhodopsin, from 0.5 kcal/mol in the A helix to 4 kcal/mol in the ED helix. Reconstructing the energy landscape as a function of molecular extension highlights that some structural elements of bacteriorhodopsin have a higher associated unfolding energy cost and improves scientific understanding of membrane protein stability. 

INTRO

(Protein folding, emphasis on landscapes and transiently occupied states)

(SMFS to estimate landscapes, higher time resolution)

(Landscape reconstruction techniques, full versus single-barrier)

METHODS

(Discuss science paper, Figure 1)

(Discuss full energy landscape reconstruction)


----
![](./Figures/combined.png)
{#label_fig:diagram} **Bacteriorhodopsin's energy landscape determined by single-molecule force spectroscopy.** (**A**) A cartoon of the pulling geometry used to unfold bR. (**B**) A representative force versus molecular extension curve, or force-extension curve. The first 20nm of all curves are not analyzed due to non-specific surface adhesion. The colored bars correspond to the regions where the helical pairs from (A) rupture. Many such curves are used for reconstruction an energy landscape via an inverse Weierstrass transform. (**C**) A detailed plot of the A helix, which exhibits transitions between states. (**D**) A plot detailing many force extension curves of bR with worm-like chain polymer fits overlayed onto major states. (**E**) An equilibrium assay demonstates rapid transitions between intermediate states in the E helix. (**F**) The energy landscape obtained via applying  $p_{\text{fold}}$ to equilibrium data as in (E) (XXX). Figures (D-F) are reproduced with permission from XXX.
----

Jarzynski's equality (@Jarzynski1996Nonequilibrium) is a thermodynamic relationship between the Helmholtz free energy of a system, $A(z)$, along a reaction coordinate $z$ to the measured work done along that coordinate, $W(z)$:

$e^{-\beta A(z)} = <e^{-\beta W(z)}>$,

where the average is taken over many independent experiments, each starting and ending at the same choice of z. For AFM-SMFS, the Helmholtz energy includes the energy stored in the cantilever used to apply forces. Jarzynski's equality is exact only in the limit as $N \rightarrow \infty$, but is approximately true for finite $N$. The equality is remarkable because it relates the work done during many repetitions of a non-equilibrium process to the equilibrium free energy difference as a function of the reaction coordinate. In practice, SMFS experiments apply Jarzynski's equality by repeatedly folding or unfolding a single molecule of interest using a force probe. In this case, and the work is the integral of the force as a function of the reaction coordinate $z$, where $z$ is the position of the cantilever. For example, a linear ramp moving from $z_0$ with constant velocity $v_0$ would have $z(t)=v_0 t + z_0$.

The effect of the force probe on the energy lanscape must be removed for an accurate determination of an energy landscape and the molecular extension. The inverse Weierstrass transform determines molecular free energy as a function of molecular extension by correcting for the perturbation of the force probe. The correction modifies Jarzynski's equality assuming a stiff, harmonic pulling apparatus (@Hummer2010_free). The weighted histogram analysis method (@Minh2008_optimized) can be used when the probe is not a harmonic spring (*e.g.* with DNA molecules linking the probe to the system of interest) or when the probe stiffness, including possible linkers, is not much greater than the system stiffness.

----
![](./Figures/iwt_diagram.png)
{#label_fig:full} **Energy landscape reconstruction of bacteriorhodopsin reveals significant intra-molecule variation in unfolding energy.**  (**A**) A heat map of all force-extension curves used in this work. Data within 20nm of the surface are excluded due to surface adhesion. (**B**) The Helmholtz free energy, the red, dashed and dotted line $A(z)$, was corrected as in @Hummer2010_free to obtain the free energy as a function of molecular extension, the black dotted line $\Delta G$. The derivative and second derivative of the Helmholtz free energy, denoted by a dot, are with respect to extension. (**C**) The $\Delta G$ from (B), where the shaded region gives the standard deviation from three non-overlapping subsets of the data in (A). 
----

Figure 2 demonstrates the inverse Weierstrass as applied to bacteriorhodopsin. Figure 2a is a heatmap of 168 force-extension curves. This heatmap represents the ensemble of measurements needed for applying Jarzynski's inequality and the corrections of the inverse Weierstrass transform. Figure 2b shows A(z) as obtained by Jarzynski's inequality (@Jarzynski1996Nonequilibrium), as well as the inverse Weierstrass corrections which are functions of A(z), the stiffness $k$, and the temperature. The corrections were obtained as in @Hummer2010_free and @Mihn2008_Free. 

DISCUSSION / RESULTS


----
![](./Figures/gallery.png)
{#label_fig:helixE} **Unfolding energies vary within a single bacteriorhodopsin molecule.** This figure shows the energy landscape applied to each helical region of bacteiorhodopsin separately. The left axis is the change in free energy, plotted as a dotted line which is color-coded by the helix. The standard deviation is shown as a shaded region around the mean. Errors were determined by applying the inverse Weierstrass transform to each of three equal-sized subsets. The right axis denotes free energy change, denoted by a thick purple line. 
----


(Motivation for full energy landscape reconstruction) 

(Discuss full energy landscape of bR -- compare average energy per amino acid and full landscape to previous results. )


The local stiffness of bacteriorhodopsin introduced error in the reconstruction of the ED helix. The inverse Weierstrass requires that stiffness of the probe be much greater than the stiffness of the reconstructed landscape. As shown in {#ref_fig:helixE}, the greatest stiffness of the bacteriorhodopsin, at the top of the ED helix, is about 8 kcal/(mol$\cdot$ $\text{nm}^2$)$\approx$ 50pN/nm, larger than the cantilever stiffness of 20pN/nm. Therefore, the top of the ED Helix was likely poorly reconstructed by the inverse Weierstrass and represents a lower bound on the true landscape. For the CB and A helices of the protein, where the stiffness was at least an order of magnitude lower, the reconstruction of the landscape was better (XXX cite). In addition, the almost-neglible correction to the landscape from the $\ddot{A}$ term in the inverse Weierstrass transform outside of the ED Helix (Figure {#label_fig:full}) confirmed that higher-order corrections were unlikely to effect the landscape of the CB and A helices.

(Discuss using greater stiffnesses) 

(Discuss per-helix results)

CONCLUSIONS 

