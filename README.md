# Installment 

1. edit 'src/Makefile'

2. set up fortran90 compiler (F90): ifort, gfortran, or others

3. set up the library of LAPACK:  MKL, ACML, ATLAS, or lapack

4. make

5. after compiling successfully, the following executable files are generated: bs_kpoints_gen, bandstructureplot, fkpt_vbm_cbm,  em_kpoints_gen, em_linefit, fs_bxsf

The usage of the above executable files is listed in README.usage.

# Note
- Recommend the output files of the self-consistent field (SCF) calculations using the tag 'ISMEAR=-5' to be saved in a directory named 'scf/'. 

- Recommend the output files of band-structure calculations to be saved separately in a directory named 'band/'. So when one runs the executable file 'bandstructureplot' in the direcctory of 'band/', the value of Fermi level will be automatically extracted from the 'OUTCAR' file in the '../scf/' directory.
