
**Free Energy Landscape of Bacteriorhodopsin Reveals Local Variations in Unfolding Energy**

P. Heenan, H. Yu, M. Siewny , T. Perkins

# ABSTRACT 

Precisely quantifying the energetics the drive the folding of membrane proteins into a lipid bilayer remains challenging. More than 15 years ago, atomic force microscopy (AFM) emerged as a powerful tool to mechanically extract individual proteins from a lipid bilayer. Concurrently, fluctuation theorems, such as the Jarzynski equality, were applied to deduce equilibrium free-energies ($\Delta G_0$) from non-equilibrium single-molecule force spectroscopy (SMFS) records. These two advances in single-molecule studies determined the free-energy of the model membrane protein bacteriorhodopsin in its native lipid bilayer. To elucidate this free-energy landscape at a higher resolution, we applied two recent developments. First, as an input to the reconstruction, we acquired force-extension curves with a 100-fold higher time resolution and 10-fold higher force precision than traditional AFM studies of membrane proteins. Next, by using the inverse Weierstrass transform of the Jarzynski free-energy integral, we avoided convolving the free energy associated with the force probe with the landscape of the molecule under study, bacteriorhodopsin. The resulting landscape with a ~2 Å resolution yielded an average unfolding free-energy per amino acid (aa) of 0.7 $\pm$ 0.1 kcal/mol, in agreement with past bulk and single-molecule studies. This high-resolution landscape derived from a non-equilibrium measurement also agreed with the one prior equilibrium measurement of $\Delta G_0$ for a membrane protein, which investigated a particular three aa transition and yielded 2.7 kcal/mol/aa.  Hence, while average unfolding $\Delta G_0$ per aa is a useful metric, our high-resolution study highlighted significant local variation from the mean.

# INTRO

Determining the structural causes of membrane protein stability is critical due to the common role membrane proteins play as drug targets.  In cells, membrane proteins are embedded in a lipid bilayer, and quantifying the intra- and inter-molecular forces driving membrane protein folding into a bilayer remains challenging. Single molecule techniques have emerged as powerful tools for studying membrane proteins in their native lipid bilayers and determining molecular stability and folding pathways by energy landscape reconstruction [@stanley_process_2008].

The energy landscape is a powerful theoretical framework for understanding protein folding and unfolding.  Structural states are represented by minima in an energy landscape, and protein folding occurs as a path along the landscape. Although many paths are possible due to local minima in the landscape and the stochastic nature of protein folding, not all pathways are equally likely, and molecules may spend only microseconds in an unstable state before transitioning. Measuring these short-lived intermediate states along the folding pathway is needed for detailed insights into the energy landscape. In order to determine properties of the free energy landscape, single molecule force spectroscopy (SMFS) techniques repeatedly unfold membrane proteins. In atomic force microscopy (AFM) unfolding experiments, a cantilever with a nanometer-scale tip attaches to and dissociates a surface-bound molecule. AFM experiments have measured the folding pathways of the membrane protein bacteriorhodopsin by measuring barriers between states, distances between states, and even the protein's entire molecular energy landscape [@janovjak_valleys_2008].  In spite of these advances, elucidating fine-grained details of the landscape requires measuring 

Using AFM to measure rarely occupied states during an unfolding pathway has traditionally been problematic due to poor time resolution [@janovjak_valleys_2008]. Although energy barriers between major states of bacteriorhodopsin have been measured using AFM as a function of point mutation [@janovjak_valleys_2008]; C- or N-terminal pulling geometry [@kessler_unfolding_2006]; and environmental conditions [@park_stabilizing_2007;preiner_free_2007], these studies either failed to reconstruct the entire energy landscape; were limited by poor temporal resolution; or failed to remove the perturbing effect of the force probe on the energy landscape. However, quantifying the location and lifetime of an intermediates as a function of reaction coordinate is required for a detailed understanding of the energy landscape's shape and roughness.

Identification of intermediate protein unfolding states has been accelerated by order-of-magnitude advances in AFM cantilever force precision and time resolution [@yu_hidden_2017]. Ten-fold higher force precision and 100-fold improvements in time resolution are facilitated by removing the gold coating of certain commercially available cantilevers and using a focused ion-beam to modify the cantilever geometry [@sullan_atomic_2013;edwards_optimizing_2017]. By taking advantage of these improvements, our experiments reconstructed the unfolding energy landscape of bacteriorhodopsin in its native lipid bilayer as a function of molecular extension by repeated mechanical dissociation (see Figure {#ref_fig:diagram}a). The reconstruction applied the inverse Weierstrass transform of the Jarzynski free-energy integral [@jarzynski_nonequilibrium_1997], which relates the free energy landscape to an ensemble of non-equilibrium SMFS experiments.  This work reports the first energy landscape of bacteriorhodopsin which is reconstructed with single-amino acid resolution and which removes the effect of the force probe on the energy landscape. 


# METHODS

(Discuss science paper, Figure 1)

## Experimental

The details surrounding cantilever modification and sample preparation are discussed elsewhere [@yu_hidden_2017]. Briefly...

(Discuss full energy landscape reconstruction)

----
![](./Figures/combined.png)
{#label_fig:diagram} **Bacteriorhodopsin's energy landscape determined by single-molecule force spectroscopy.** (**A**) A cartoon of the pulling geometry used to unfold bR. (**B**) A representative force versus molecular extension curve, or force-extension curve. The first 20nm of all curves are not analyzed due to non-specific surface adhesion. The colored bars correspond to the regions where the helical pairs from (A) rupture. Many such curves are used for reconstruction an energy landscape via an inverse Weierstrass transform. (**C**) A detailed plot of the A helix, which exhibits transitions between states. (**D**) A plot detailing many force extension curves of bR with worm-like chain polymer fits overlayed onto major states. (**E**) An equilibrium assay demonstrates rapid transitions between intermediate states in the E helix. (**F**) The energy landscape obtained via applying  $p_{\text{fold}}$ to equilibrium data as in (E) (XXX). Figures (D-F) are reproduced with permission from XXX.
----

## Landscape reconstruction

Jarzynski's equality (@Jarzynski1996Nonequilibrium) is a thermodynamic relationship between the Helmholtz free energy of a system, $A(z)$, along a reaction coordinate $z$ to the measured work done along that coordinate, $W(z)$:

$e^{-\beta A(z)} = <e^{-\beta W(z)}>$,

where the average is taken over many independent experiments, each starting and ending at the same choice of z. For AFM-SMFS, the Helmholtz energy includes the energy stored in the cantilever used to apply forces. Jarzynski's equality is exact only in the limit as $N \rightarrow \infty$, but is approximately true for finite $N$. The equality is remarkable because it relates the work done during many repetitions of a non-equilibrium process to the equilibrium free energy difference as a function of the reaction coordinate. In practice, SMFS experiments apply Jarzynski's equality by repeatedly folding or unfolding a single molecule of interest using a force probe. In this case, and the work is the integral of the force as a function of the reaction coordinate $z$, where $z$ is the position of the cantilever's base. For example, a linear ramp moving from $z_0$ with constant velocity $v_0$ would have $z(t)=v_0 t + z_0$.

The effect of the force probe on the energy landscape must be removed for an accurate determination of an energy landscape and the molecular extension. The inverse Weierstrass transform of Jarzynski's free-energy integral, hereafter referred to as the inverse Weierstass transform, determines molecular free energy as a function of molecular extension by correcting for the perturbation of the force probe. The correction modifies Jarzynski's equality assuming a stiff, harmonic pulling apparatus (@Hummer2010_free). The weighted histogram analysis method (@Minh2008_optimized) can be used when the probe is not a harmonic spring (*e.g.* with DNA molecules linking the probe to the system of interest) or when the probe stiffness, including possible linkers, is not much greater than the system stiffness.

----
![](./Figures/iwt_diagram.png)
{#label_fig:full} **Energy landscape reconstruction of bacteriorhodopsin reveals significant intra-molecule variation in unfolding energy.**  (**A**) A heat map of all force-extension curves used in this work. Data within 20nm of the surface are excluded due to surface adhesion. (**B**) The Helmholtz free energy, the red, dashed and dotted line $A(z)$, was corrected as in @Hummer2010_free to obtain the free energy as a function of molecular extension, the black dotted line $\Delta G$. The derivative and second derivative of the Helmholtz free energy, denoted by a dot, are with respect to extension. (**C**) The $\Delta G$ from (B), where the shaded region gives the standard deviation from three non-overlapping subsets of the data in (A). 
----

Figure 2 demonstrates the inverse Weierstrass as applied to bacteriorhodopsin. Figure 2a is a heatmap of 168 force-extension curves. This heatmap represents the ensemble of measurements needed for applying Jarzynski's inequality and the corrections of the inverse Weierstrass transform. Figure 2b shows A(z) as obtained by Jarzynski's inequality (@Jarzynski1996Nonequilibrium), as well as the inverse Weierstrass corrections which are functions of A(z), the stiffness $k$, and the temperature. The corrections were obtained as in @Hummer2010_free and @Mihn2008_Free. 

# DISCUSSION / RESULTS


----
![](./Figures/gallery.png)
{#label_fig:helixE} **Unfolding energies vary within a single bacteriorhodopsin molecule.** This figure shows the energy landscape applied to each helical region of Bacteriorhodopsin separately. The left axis is the change in free energy, plotted as a dotted line which is color-coded by the helix. The standard deviation is shown as a shaded region around the mean. Errors were determined by applying the inverse Weierstrass transform to each of three equal-sized subsets. The right axis denotes free energy change, denoted by a thick purple line. 
----


(Motivation for full energy landscape reconstruction) 

(Discuss full energy landscape of bR -- compare average energy per amino acid and full landscape to previous results. )

The unfolding energy of bacteriorhodopsin depends on the structural element being unfolded. Figure {#ref_fig:helixE} reports the energy landscape and unfolding energy per amino acid as a function of molecular extension for the ED, CB, and A helices. Of these, the ED helix has the highest associated unfolding energy per amino acid, followed by the CB and A helices. Although this energy landscape provides a good estimation for the CB and A helices, the reconstruction of the ED helical pair is limited by its extremely high stiffness. 

The large stiffness of the ED helical pair of bacteriorhodopsin, compared to the force probe, introduced error in its energy landscape reconstruction. The inverse Weierstrass requires that stiffness of the probe be much greater than the stiffness of the reconstructed landscape. As shown in {#ref_fig:helixE}, the greatest stiffness of the bacteriorhodopsin, at the top of the ED helix, is about 8 kcal/(mol$\cdot$ $\text{nm}^2$)$\approx$ 50pN/nm, larger than the cantilever stiffness of 20pN/nm. Therefore, the top of the ED Helix was likely poorly reconstructed by the inverse Weierstrass and represents a lower bound on the true landscape. For the CB and A helices of the protein, where the stiffness was at least an order of magnitude lower, the reconstruction of the landscape was better (XXX cite). In addition, the almost-negligible correction to the landscape from the $\ddot{A}$ term in the inverse Weierstrass transform outside of the ED Helix (Figure {#ref_fig:full}) confirmed that higher-order corrections were unlikely to effect the landscape of the CB and A helices.

# CONCLUSION

[xxx bridge...]

The reconstruction presented provides insight into the unfolding energies ... . The


(Discuss using greater stiffnesses) 

(Discuss per-helix results)

(Discuss applying to environmental, other puling methods, site-specific (GF helix), retinol, etc etc.)

Further work: 

* unfolding,refolding.
* dependence on conditions 
* comparison with equilibrium measurements
* removing the effect of handles 

