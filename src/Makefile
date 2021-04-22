#note the option for the fort compiler
F90=ifort
#F90=gfortran

#LAPACKLIBS=-llapack -lblas #  [for Mac os] comment this line if you are not using LAPACK.
                             # On systems with a FORTRAN (not f2c'ed) version of LAPACK, -lf2c is
                             # not necessary; on others, -lf2c is equivalent to -lF77 -lI77

#LAPACKLIBS=-L/usr/local/atlas/lib -llapack -lcblas -lf77blas -latlas -lf2c # This works with   the ATLAS updated lapack and Linux_P4SSE2
                                    # from   http://www.netlib.org/atlas/archives/linux/
#LAPACKLIBS=-L/opt/local/lib -llapack -lblas -latlas


#Use acml
#LAPACKLIBS=/home/z-hou/soft/acml531/gfortran64/lib/libacml.a
#LAPACKLIBS=/work/z-hou/soft/acml531/ifort64/lib/libacml.a
#LAPACKLIBS=-L/admin/app/gnu/lapack-3.5.0/lib -llapack -lblas 
#LAPACKLIBS=-L/opt/local/lib -llapack -lblas 

#Use MKL
#MKLROOT =/admin/app/intel/mkl
MKLROOT =/export/apps/intel/composer_xe_2013_sp1.0.080/mkl

MKL_PATH=$(MKLROOT)/lib/intel64
LAPACKLIBS=-L$(MKL_PATH)  -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread
#or
#LAPACKLIBS=-L/opt/intel/Compiler/11.1/046/mkl/lib/em64t -lmkl_lapack -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread

LIBS=$(LAPACKLIBS)

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<
	
#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------
SRC = precmod.f90 constant.f90 stringmod.f90 readinpem.f90 readposcar.f90  \
      eigenvalueparam.f90   readprocar.f90 getefermi.f90\
      bandstructureplot.f90  highsymkpoints.f90  writehighk.f90 \
      writegnuplotemplate.f90 getkdistances.f90 bs_kpoints_gen.f90 \
      fitting.f90 findinv.f90 em_linefit.f90  em_kpoints_gen.f90  \
      fkpt_vbm_cbm.f90 getlatinfo.f90 fs_bxsf.f90 

OBJ = $(SRC:.f90=.o)

all:  $(OBJ)   em_linefit em_kpoints_gen bandstructureplot  bs_kpoints_gen fkpt_vbm_cbm   getlatinfo fs_bxsf
	

bandstructureplot:	$(OBJ)
	$(F90) $(F90_OPTS) -o bandstructureplot bandstructureplot.o  precmod.o constant.o getkdistances.o findinv.o readprocar.o eigenvalueparam.o stringmod.o  readposcar.o getefermi.o

fs_bxsf:	$(OBJ)
	$(F90) $(F90_OPTS) -o fs_bxsf fs_bxsf.o  precmod.o constant.o findinv.o eigenvalueparam.o stringmod.o  readposcar.o getefermi.o

bs_kpoints_gen:	$(OBJ)
	$(F90) $(F90_OPTS) -o bs_kpoints_gen bs_kpoints_gen.o  precmod.o constant.o getkdistances.o findinv.o readposcar.o   highsymkpoints.o  writehighk.o writegnuplotemplate.o

getlatinfo:	$(OBJ)
	$(F90) $(F90_OPTS) -o getlatinfo getlatinfo.o  precmod.o constant.o  findinv.o readposcar.o

fkpt_vbm_cbm:	$(OBJ)
	$(F90) $(F90_OPTS) -o fkpt_vbm_cbm fkpt_vbm_cbm.o  precmod.o constant.o   eigenvalueparam.o stringmod.o 

em_kpoints_gen:	$(OBJ)
	$(F90) $(F90_OPTS) -o em_kpoints_gen em_kpoints_gen.o findinv.o readposcar.o readinpem.o precmod.o constant.o

em_linefit: $(OBJ)  
	$(F90)  $(F90_OPTS)  -o em_linefit em_linefit.o  precmod.o constant.o fitting.o  findinv.o  readposcar.o readinpem.o eigenvalueparam.o  stringmod.o  $(LIBS) 

cleanall:
	@rm -f *.o *.mod em_kpoints_gen  em_linefit  bandstructureplot  bs_kpoints_gen fkpt_vbm_cbm   getlatinfo
clean:
	@rm -f *.o *.mod 

#end
