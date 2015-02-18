# This is a VERY striped down makefile for compiling sac with gfortran

# Check to see if the ENV variables are defined for the compiler
ifndef F90
    F90=mpif90
endif
ifndef F90FLG
    F90FLG=-std=f2008
endif
F90FLAGS= $(F90FLG)

FOR=$(F90)
FORFLG=$(F90FLAGS)
PREFOR=$(F90TOF90)

VACPP_=./vacpp.pl
VACPP =./vacpp.pl

BLASLIB = blas
LAPACKLIB = lapack

F=.f90
O=.o

VACDIR=.

################# Definitions for source files #############################

LIBS = vacdef$F
INCLUDES = vacdef$F vacpar$F vacusrpar$F

VAC_FOR = vac$F vacio$F vacgrid$F vacphys0$F vacphys$F vacusr$F
VAC_OBJ = vac$O vacio$O vacgrid$O vacphys0$O vacphys$O vacusr$O

# The VACCD, VACMC, VACFCT, VACTVDLF, VACTVD, VACIMPL, VACPOISSON, VACMPI 
# modules are removed or added in the following two lines by "setvac". 
# DO NOT TOUCH THESE TWO LINES:
VACFOR = $(VAC_FOR) vaccd$F
VACOBJ = $(VAC_OBJ) vaccd$O

ROEOBJ = roetest$O vacphys0$O vacphys$O

################### Translation and dependency rules #######################

PREPROC= $(VACPP_) 

.SUFFIXES:
.SUFFIXES: .t $F  $O

$(VACFOR) vacini$F vacdef$F vacpar$F vacusrpar$F        : $(PREPROC)
roetest$F                           : $(PREPROC)

$(VACOBJ) vacini$O vacall$O vaciniall$O vacsmall$O roetest$O    : $(INCLUDES)

vacusrpar$F: vacusrpar.t
	$(VACPP) $< $(PREFOR) > $@
	$(FOR) -c vacdef.f90

vacpar$F: vacpar.t
	$(VACPP) $< $(PREFOR) > $@

vacdef$F: vacdef.t
	$(VACPP) $< $(PREFOR) > $@


# General precompilation rule
.t$(F):
	$(VACPP) $< $(PREFOR) > $@

# General Fortran compilation rule
$(F)$(O):
	$(FOR) $(FORFLG) -c $<



########### Extra dependencies for some files with "INCLUDE:" statements
vacphys$F: vacphys.t vacphys.mhd0.t vacphys.mhdres.t
	$(VACPP) $< $(PREFOR) > $@

vacgrid$F: vacgrid.t
	$(VACPP) $< $(PREFOR) > $@

vac :$(VACOBJ) $(LIBS_)
	$(FOR) $(FORFLG) -o $(VACDIR)/vac $(VACOBJ) $(LIBS)

###### Removing object files, precompiled Fortran files, and symbolic links

clean   :
	rm -f vac*$F roetest$F Vvac*$F vac*.d *$O *~ *mod

cleanall: clean
	rm -f vacphys.t vacphys0.t vacpar.t vacusrpar.t vacusr.t 
