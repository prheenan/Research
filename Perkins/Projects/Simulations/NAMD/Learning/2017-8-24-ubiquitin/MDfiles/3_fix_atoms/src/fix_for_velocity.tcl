# script to automate setting up molecule for fixed velocity pulling

mol new \
    ../1_solvation/ubq.psf
mol addfile \
    ../2_equilibration/src/ubq_ww_eq.pdb

set allatoms [atomselect top all]
# make one atom fixed
$allatoms set beta 0
set fixedatom [atomselect top "resid 1 and name CA"]
# keep the first residue alpha carbon to be fixed 
$fixedatom set beta 1
# mkae an atom to pull on (last residue alpha carbon)
$allatoms set occupancy 0
set smdatom [atomselect top "resid 76 and name CA"]
$smdatom set occupancy 1
$allatoms writepdb fixed_ubq.ref

set smdpos [lindex [$smdatom get {x y z}] 0]
set fixedpos [lindex [$fixedatom get {x y z}] 0]
vecnorm [vecsub $smdpos $fixedpos]
