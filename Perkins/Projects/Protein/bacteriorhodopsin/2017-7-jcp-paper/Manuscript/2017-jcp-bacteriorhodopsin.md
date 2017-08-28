**Reconstructing Order-of-Magnitude Unfolding Energy Variations in Bacteriorhodopsin's Residues by Single Molecule Force Spectroscopy**

P. Heenan, H. Yu, M. Siewny , T. Perkins


INTRO

(Protein folding, emphasis on landscapes and transiently occupied states)

(SMFS to estimate landscapes, higher time resolution)

(Landscape reconstruction techniques, full versus single-barrier)

METHODS

(Discuss science paper, Figure 1)

(Discuss full energy landscape resonstuction)


----
![](./Figures/diagram.png)
{#label_fig:diagram} Atomic force microscopy measures the unfolding and refolding of bacteriorhodopsin (bR) under mechanical load. **(A)** A cartoon of the pulling geometry used to unfold bR. **(B)** A representative force versus extension molecular extension curve. The first 20nm of all curves are not analyzed due to non-specific surface adhesion. The colored bars correspond to the regions where the helical pairs from (A) rupture. **(C)** A detailed plot of the A helix, which exhibits transitions between states. **(D)** A plot detailing many force extension curves of bR with worm-like chain polymer fits overlayed onto major states. **(E)** An equilibrium assay demonstates rapid transitions between intermediate states in the E helix. **(F)** The energy landscape obtained via an inverse boltzmann transform (XXX) . Figures (D-F) are reproduced with permission from XXX.
----

Jarzynski's equality (@Jarzynski1996Nonequilibrium) is a thermodynamic relationship between the Helmholtz free energy difference, $A(z)$, along a reaction coordinate $z$ to the measured work done along that coordinate, $W(z)$:

$e^{-\beta A(z)} = <e^{-\beta W(z)}>$,

where the average is taken over many independent experiments, each starting and ending at the same choice of z. Jarzynski's equality is exact only in the limit as $N \rightarrow \infty$, but is approximately true for finite $N$. The equality is remarkable because it relates the work done during many repetitions of a non-equilibrium process to the equilibrium free energy difference as a function of the reaction coordinate. In practice, SMFS experiments apply Jarzynski's equality by repeatedly folding or unfolding a single molecule of interest using a force probe. In this case, the reaction coordinate, $z$, is the end-to-end molecular distance, referred to as the extension, and the work is the integral of the force as a function of the extension.

The effect of the force probe on the energy lanscape must be removed for an accurate determination of an energy landscape. The inverse Weierstrass transform determines molecular free energy as a function of molecular extension correcting for the perturbation of a force probe. The correction modifies Jarzynski's equality assuming a stiff, harmonic pulling apparatus (@Hummer2010_free). The weighted histogram analysis method (@Minh2008_optimized) can be used when the probe is not a harmonic spring (e.g. with DNA molecules linking the probe to the system of interest) or when the probe stiffness, including possible linkers, is not much greater than the system stiffness. (XXX run this, discuss).


----
![](./Figures/iwt_diagram.png)
{#label_fig:full} Energy landscape reconstruction of bacteriorhodopsin reveals significant intra-molecule variation in unfolding energy.  **(A)** A heat map of all force-extension curves used in this work. Data within 20nm of the surface are excluded due to surface adhesion. **(B)** The Helmholtz free energy, $A(z)$, is corrected as in @Hummer2010_free to obtain the free energy as a function of molecular extension, $\Delta G$. **(C)** The mean free energy at zero force (black dotted line) and free energy change per amino acid (purple line) reconstructed using an inverse Weierstrass transform applied to the data in (A).  The shaded region gives the standard deviation from three non-overlapping subsets of the data in (A). 
----

Figure 2 demonstrates the inverse Weierstrass as applied to bacteriorhodopsin. Figure 2a is a heatmap of 168 force-extension curves. This heatmap represents the ensemble of measurements needed for applying Jarzynski's inequality and the corrections of the inverse Weierstrass transform. Figure 2b shows A(z) as obtained by Jarzynski's inequality (@Jarzynski1996Nonequilibrium), as well as the inverse Weierstrass corrections which are functions of A(z), the stiffness $k$, and the temperature. The corrections are obtained by the method of @Hummer2010_free and @Mihn2008_Free.

DISCUSSION / RESULTS


----
![](./Figures/gallery.png)
{#label_fig:helixE} The ED helical pair has a significantly higher average unfolding energy than the CB helical pair or the A Helix. This figure is formatted as {#ref_fig:full}, but applied to the region of the force-extension curve of the ED Helical pair.
----



(Motivation for full energy landscape reconstruction) 

(Discuss full energy landscape of bR -- compare average energy per amino acid and full landscape to previous results. )

(Discuss per-helix results)

CONCLUSIONS 

