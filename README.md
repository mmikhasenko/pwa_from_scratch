# PWA from Scratch

## Running analysis

1. Precalculation: `analysis/02_precalculation.jl`
2. Compute integrals: `analysis/03_integrals.jl`
3. Compute integrals: `analysis/04_integrals_fu.jl`
4. Perform fit: `analysis/05_llhfit.jl`


## Validation of integrals

If you want to compare the different methods to calculate the phase-space integrals do the following steps:

1. Get files to compare to `/mnt/data/compass/2008/integrals_Florian/` on bonn local computers

2. extract the values from it.
   For this use the c++-code "save_memory.cc".
   Inside you have to specify the indices of the offdiagonal. If you don't want to extract them, set both to the same number.
   The indices of the different waves can be found in one of the integral files.

3. run `tests/integral_validation_test.jl` in julia,
   specifying the wave indices that you want to calculate the integrals for. Here the indices are different from them in step 1. Find them out by looking at the file "wavelist_formatted.txt".
   There are different functions for diagonals and interferences (offdiagonals).

4. run `tests/integral_validation_plot.jl` in julia, specifying again the same indices as in steps 1&2.
