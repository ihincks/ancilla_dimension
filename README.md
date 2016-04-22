# Description

This set of MATLAB files was used to create the computed data appearing in the paper *“Ancilla dimension in quantum channel discrimination”* by Daniel Puzzuoli and John Watrous.

The scripts `compute_1_norms` and `compute_1_norms_herm` are the functions called to produce the data, with `InducedSchattenNormMod` and `InducedSchattenNormModHerm` being where the main functionality for computing the norms lies.

The code relies on [QETLAB](http://qetlab.com), and was written using version 0.9. It directly depends on the QETLAB functions `PartialTrace` and `PartialTranspose`, and the functions `InducedSchattenNormMod` and `InducedSchattenNormModHerm` are modified versions of the QETLAB function `InducedSchattenNorm`. See the documentation in `InducedSchattenNormMod` for details on the modifications.

The data that appears in the paper is stored in the files `norm_bounds.mat` and `norm_bounds_herm.mat`. The data is saved as a struct, and the data structure is explained in the documentation at the beginning of the script `compute_1_norms`. 
