
Scripts to obtain QM/Amber surface for methanol group rotation
in hmC. Outputs are Gaussian09 input files. 

The FF data should be from a MD simulations where the dihedral potential for
C4-C5-C55-OH5 and C5-C55-OH5-HO5 where set to zero.

do_parse_g09_plumed_plot_maps_and_projections.py reads all the g09 outputs and 
free energies from Free energy surface (typically from metadynamics), and computes
the 2D free energy surfaces. 


