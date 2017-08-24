
mol new ../1_solvation/ubq.psf
mol addfile ./src/ubq_gbis_eq.restart.coor
set selprotein [atomselect top protein]
$selprotein writepdb ./src/ubq_ww_eq.pdb
