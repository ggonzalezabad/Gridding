###########################################################
#                                                         #
#   Makefile for the SAO PGE L2 Access Software           #
#                                                         #
#   Author:  Thomas P. Kurosu                             #
#                                                         #
#            Harvard-Smithsonian Center for Astrophysics  #
#            Atomic and Molecular Physics Division        #
#            60 Garden Street (MS 50)                     #
#            Cambridge, MA 02138, USA                     #
#            Tel/Fax: +1 -- 617 - 495 7213                #
#            EMail tkurosu@cfa.harvard.edu                #
#                                                         #
#                                                         #
#   Last modified: 8 February 2005                        #
#                                                         #
###########################################################

# ================================================
# Define some environment variables that
# are unlikely to change between host-types
# ================================================
AR      = ar
ARFLAGS = ruv
MV      = mv
RM      = rm -f
SHELL   = /bin/sh


# ========================================================
# Location of HDF et al. libraries. Depends on whether we
# are compiling in 32bit or 64bit.
# ========================================================
ifeq (32,$(bitversion))
   OMIUTIL = /data/tempo2/ggonzale/TEMPO_code/ots/gnu
   FLAVOR  = linux32
else
   OMIUTIL = /data/tempo2/ggonzale/TEMPO_code/ots/gnu
   FLAVOR  = linux64
endif

# ==============
# Target Library
# ==============
TARGET = SAO_OMIL2_ReadLib

FC        = gfortran
CC        = gcc
FFLAGSopt = -Ofast -mcmodel=medium 
FFLAGSdbg = -g -C -fbacktrace -Wall -fcheck=all
CFLAGS    = -Ofast -mcmodel=medium
LNFLAGS   = -static


# ------------------------------------------------
# *.mod files are found with "-I"
# ------------------------------------------------
MODFLG     = -I

# ----------------------------- 
# Do we want DEBUG information?
# ----------------------------- 
ifeq (yes,$(debug))
   FFLAGS      = $(FFLAGSdbg)
   SLATECFLAGS = -Bstatic
else
   FFLAGS      = $(FFLAGSopt)
   SLATECFLAGS = $(FFLAGS)
endif

# ----------------------------
# Generic additional libraries
# ----------------------------
ADDLIBS     = -lz -lsz -lm -ldl
PSPLINELIBS = -lpspline -lezcdf
NETCDFLIBS  = -lnetcdf
HE4LIBS     = -lhdfeos -lmfhdf -lGctp -ldf -ljpeg
HE5LIBS     = -lhe5_hdfeos -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran -lhdf5

# =================================
# HDF5 and HDFEOS5 libraries to use
# =================================
HDFEOS_INC  = $(OMIUTIL)/hdfeos/include
HDFEOS_LIB  = $(OMIUTIL)/hdfeos/lib/$(FLAVOR)
HDFEOS5_INC = $(OMIUTIL)/hdfeos5/include
HDFEOS5_LIB = $(OMIUTIL)/hdfeos5/lib/$(FLAVOR)

IFLAGS = -I$(OMIUTIL)/include -I$(OMIUTIL)/mod \
	 -I$(HDFEOS_INC) -I$(HDFEOS5_INC) 

LFLAGS = -L$(OMIUTIL)/lib -L$(HDFEOS_LIB) -L$(HDFEOS5_LIB)

LIBS   = $(PSPLINELIBS) $(NETCDFLIBS)                \
         $(HE5LIBS) $(HE4LIBS) $(ADDLIBS)


# ================================================
# All sources are "sourced out" to an include file
# ================================================
BASMOD = SAO_OMIL2_ReadLib_Basic_module
HE5MOD = SAO_OMIL2_ReadLib_he5_module
SPHMOD = SAO_OMIL2_ReadLib_SphGeo_module
TSLMOD = SAO_OMIL2_ReadLib_Tessel_module
MISMOD = SAO_OMIL2_ReadLib_Misc_module
MEDMOD = SAO_OMIL2_median_module
NLSMOD = SAO_elsunc_fitting_module
A2HMOD = SAO_OMIL2_CreateHDF5_module
ALLMOD = SAO_OMIL2_allocation_module

#G2RMOD = SAO_GOME2L2_ReadLib_he5_module

CSOURCES     = unbufferSTDOUT.c

F90SOURCES = $(BASMOD).f90     \
             $(ALLMOD).f90     \
             $(HE5MOD).f90     \
             $(NLSMOD).f90     \
             $(SPHMOD).f90     \
             $(TSLMOD).f90     \
             $(MISMOD).f90     \
             $(MEDMOD).f90     \
             $(A2HMOD).f90     \
             ezspline_interpolation.f90 \
             SLATEC_davint.f90
#             $(G2RMOD).f90     \

OBJECTS = $(CSOURCES:.c=.o) $(F90SOURCES:.f90=.o)
OBJECTDIR = ./objects
LINKOBJECTS = $(addprefix $(OBJECTDIR)/, $(OBJECTS))

# -----------------------------------
# The PATH to search for dependencies
# -----------------------------------
VPATH = ./:$(OBJECTDIR)

# =================================================
#  Everything we want to make in this distribution
# =================================================
all : makefolders $(OBJECTS) $(target)

$(target): $(target).f90 $(OBJECTS)
	$(FC) $(target).f90 -o $(target).x $(LNFLAGS) $(FFLAGS) $(IFLAGS) $(LFLAGS) $(LINKOBJECTS) $(LIBS)

.PHONY : clean cleanobs cleanmod cleantarget cleanfolder makefolders

# ------------
# Make folders
# ------------
makefolders:
	@echo 'Creating folder object ' $(OBJECTDIR)/ '...'
	@mkdir -p $(OBJECTDIR)

# ---------------------------------------------------
#  Just in case, we want to CLEAN up once in a while
# ---------------------------------------------------
cleanobs: 
	@echo ' * Removing complied object files in ./ and' $(OBJECTDIR)/ '...'
	@-rm -f ./$(OBJECTDIR)/* ./*.o
cleanmod: 
	@echo ' * Removing compiled F90 Module files ...'
	@-rm -f *.mod *.M *genmod*
cleantarget: 
	@if [ -e ./$(target).x ]; then                             \
            echo ' * Removing compiled target' $(target).x '...';  \
            rm -f $(target).x;                                     \
         fi
cleanfolder: 
	@if [ -e ./$(OBJECTDIR) ]; then                             \
            echo ' * Removing compiled object folder' $(OBJECTDIR) '...';  \
            rm -fr $(OBJECTDIR);                                     \
         fi
clean: cleanobs cleanmod cleantarget cleanfolder

# ==================================
#  PATTERN rule for F90 compilation
# ==================================
.SUFFIXES: .f90 .f


SLATEC_davint.o:
	$(FC) -c $(SLATECFLAGS) $(IFLAGS) $<
	$(MV) $*.o $(OBJECTDIR)/

.f90.o :
	$(FC) -c $(FFLAGS) $(IFLAGS) $<
	$(MV) $*.o $(OBJECTDIR)/

.f.o :
	$(FC) -c $(FFLAGS) $(IFLAGS) $<
	$(MV) $*.o $(OBJECTDIR)/

.c.o:
	$(CC) -c $(CFLAGS) $<
	$(MV) $*.o $(OBJECTDIR)/

# ----------------------------------------------
# Short-hand notations for compiled MODULE files
# ----------------------------------------------
BASM = $(BASMOD).o
ALLM = $(ALLMOD).o
HE5M = $(HE5MOD).o
SPHM = $(SPHMOD).o
TSLM = $(TSLMOD).o
MSLM = $(MISMOD).o
MEDM = $(MEDMOD).o
NLSM = $(NLSMOD).o
A2HM = $(A2HMOD).o
G2RM = $(G2RMOD).o

# =====================================================
#  Various DEPENDENCIES, which tell Make  to recompile
#  subroutines when the source files have changed
# =====================================================
$(BASMOD).o : $(BASMOD).f90
$(ALLMOD).o : $(ALLMOD).f90   $(BASM)
$(HE5MOD).o : $(HE5MOD).f90   $(BASM)
$(SPHMOD).o : $(SPHMOD).f90   $(BASM)
$(TSLMOD).o : $(TSLMOD).f90   $(BASM)
$(MISMOD).o : $(MISMOD).f90   $(BASM)
$(MEDMOD).o : $(MEDMOD).f90   $(BASM)
$(A2HMOD).o : $(A2HMOD).f90   $(BASM)
$(NLSMOD).o : $(NLSMOD).f90   $(BASM)
$(G2RMOD).o : $(G2RMOD).f90   $(BASM)
ezspline_interpolation.o : ezspline_interpolation.f90 $(BASM)
SLATEC_davint.o          : SLATEC_davint.f90          $(BASM)

unbufferSTDOUT.o : unbufferSTDOUT.c

# =================
#  End of Makefile
# =================
