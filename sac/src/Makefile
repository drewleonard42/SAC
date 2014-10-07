# This is a VERY striped down makefile for compiling sac with gfortran
# This will become a SAC only f2py makefile hopefully 
ifndef F90
    F90=h5pfc
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
GDFDIR=./fgdfio/lib

################# Definitions for source files #############################

LIBS = vacdef$F sacgdf$F
INCLUDES = vacdef$F vacpar$F vacusrpar$F sacgdf$F
GDF_INCLUDES = $(GDFDIR)/gdf_types.F90 $(GDFDIR)/helpers_hdf5.F90 $(GDFDIR)/grid_data_format.F90 $(GDFDIR)/gdf_datasets.F90

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

$(VACFOR) vacini$F vacdef$F vacpar$F vacusrpar$F sacgdf$F       : $(PREPROC)
roetest$F                           : $(PREPROC)

$(VACOBJ) vacini$O vacall$O vaciniall$O vacsmall$O roetest$O sacgdf$0   : $(INCLUDES)

vacusrpar$F: vacusrpar.t
	$(VACPP) $< $(PREFOR) > $@
	$(FOR) -c vacdef.f90

vacpar$F: vacpar.t
	$(VACPP) $< $(PREFOR) > $@

vacdef$F: vacdef.t
	$(VACPP) $< $(PREFOR) > $@

sacgdf$F: sacgdf.t
	$(VACPP) $< $(PREFOR) > $@
	$(FOR) -c sacgdf.f90

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

gdf :
	$(FOR) $(FORFLG) -c $(GDF_INCLUDES)

vac : gdf $(VACOBJ) $(LIBS_)
	$(FOR) $(FORFLG) -o $(VACDIR)/vac $(GDF_INCLUDES) $(VACOBJ) $(LIBS) 

###### Removing object files, precompiled Fortran files, and symbolic links

clean   :
	rm -f vac*$F sac*$F roetest$F Vvac*$F vac*.d *$O *~ *.mod

cleanall: clean
	rm -f vacphys.t vacphys0.t vacpar.t vacusrpar.t vacusr.t 
