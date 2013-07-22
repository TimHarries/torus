#-------------------------------------------------------------------------
# Usage :: 
#    1. make             : Make torus with optimization flag(s)
#    2. make clean       : Clean up directory
#    3. make cfitsio=no  : Use cfitsio stubs rather than library
#    4. make debug=yes   : Make torus with debug option,
#    5. make profile=yes : Make torus with profiler option
#    6. make static=yes  : Make torus with static linking
#    7. make trace=yes   : Make torus with ITAC profiling
#    8. make openmp=yes  : Make torus with OpenMP 
#    9. make coverage=yes: Generate coverage statistics with gnu compilers 
#   10. make codecov=yes : Generate coverage statistics with ifort
#   11. make zlib=yes:   : Use z-lib compression libary to produce vtu files
#-------------------------------------------------------------------------

USEMKL  = no
cfitsio = yes
hdf = no
memcheck = yes
zlib = yes
debug   = no
profile = no
static  = no
trace   = no
openmp  = no 
v2      = yes
coverage = no
codecov  = no
withsphng = no
hydro = yes
molecular = yes
photoion = yes
sph = yes
atomic = yes
xray = yes
# These are set later on, don't set in the make command
stateq = no
datacube = no
# Run script to get subversion version number
getsvnver = yes

KNOWN_SYSTEM = no
MPIFLAG = 
#
# Set system-dependent variables
#

########################
# NAG Fortran compiler #
########################

ifeq ($(SYSTEM), nagfor)
    F95COMP  = nagfor
    F95FLAG  = -colour -kind=byte -fpp -dcfuns -f2003 -DUSEIEEEISNAN -DUSEUNIXENV
    DEBUGFLAG= -C=all -g -gline -info -nocheck_modtime
    FASTFLAG = -O4
    CFITSIOLIBS = -L/Users/${USER}/cfitsio_nagfor/lib -lcfitsio 
    ZLIBLIBS = -lz	
    EXE = torus.nagfor
    KNOWN_SYSTEM = yes
    memcheck = no
endif

####################
# g95 based builds #
####################

ifeq ($(SYSTEM), g95)
    F95COMP  = g95
    DEBUGFLAG= -g3 -fbounds-check -ftrace=full -fimplicit-none -freal=nan -Wall -Wunused-internal-procs -Wunused-parameter -Wno=140,141,112 -Werror
#    FASTFLAG = -O3 -mcpu=G5
    FASTFLAG = -O3
#    FASTFLAG = -O -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math	
    CFITSIOLIBS = -lcfitsio
    ZLIBLIBS = -lz	
    EXE = torus.g95
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), macmpi)
    F95COMP  = g95
    F95FLAG  = -fno-second-underscore
    DEBUGFLAG= -g
    STATICFLAG = -static
#    FASTFLAG = -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    MPI_INC  = -I/Users/th/macmpi
    EXE = torus.macmpi
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), ompi)
    F95COMP  = mpif90
    F95FLAG  = -Wno=155
    DEBUGFLAG= -g3 -ftrace=full -fbounds-check -fimplicit-none -freal=nan -Wall -Wunused-internal-procs -Wno=140,141,112 -Werror
    FASTFLAG = -O3 -march=pentium4 -mtune=pentium4 -malign-double
    CFITSIOLIBS = -lcfitsio
    ZLIBLIBS = -lz
    EXE = torus.ompi
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

#########################
# gfortran based builds #
#########################

ifeq ($(SYSTEM), gfortran)
    ifeq ($(mpi), yes)
        F95COMP = mpif90
        MPIFLAG = -DMPI
    else
        F95COMP = gfortran
    endif
     ZLIBLIBS = -lz
    DEBUGFLAG= -ggdb -fbacktrace -fcheck=all -Wunderflow -Wall -Werror -pedantic-errors -std=f2003 -fall-intrinsics -ffpe-trap=invalid,zero,overflow #,underflow,denormal
    FASTFLAG = -O3
#    FASTFLAG = -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    EXE = torus.gfortran
    KNOWN_SYSTEM = yes
    CFITSIOLIBS = -L/Users/${USER}/cfitsio/lib -lcfitsio
endif

# MPI build with the gfortran compiler
ifeq ($(SYSTEM), ompiosx)
    F95COMP  = mpif90 
    DEBUGFLAG= -ggdb -fbacktrace -fcheck=all -Wunderflow -Wall -Werror -pedantic-errors -std=f2003 -fall-intrinsics -ffpe-trap=invalid,zero,overflow #,underflow,denormal
    FASTFLAG = -O3
    CFITSIOLIBS = -L/Users/${USER}/cfitsio/lib -lcfitsio
    EXE = torus.ompiosx
    ZLIBLIBS = -lz	
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

######################
# ifort based builds #
######################

ifeq ($(SYSTEM), ifort)
    ifeq ($(mpi), yes)
        F95COMP = mpiifort
        MPIFLAG = -DMPI
    else
        F95COMP  = ifort 
    endif
    DEBUGFLAG= -g -O0 -warn -traceback -fpe0 -fp-stack-check -gen-interfaces -check all -check noarg_temp_created
    CFITSIOLIBS = -lcfitsio
    ZLIBLIBS = -lz
    EXE = torus.ifort
    KNOWN_SYSTEM = yes
endif

###########################
# Machine specific builds #
###########################

# The Dirac system at Leicester AKA complexity cluster
# Default is an MPI build, use mpi=no for openmp/serial
# Use 'module load cfitsio' to make FITS library available
ifeq ($(SYSTEM), complexity)
    ifeq ($(mpi), no)
       F95COMP= ifort
       MPIFLAG=
    else
       F95COMP= mpif90
       MPIFLAG= -DMPI
    endif
    FASTFLAG= -O3 
# For debugging with optimisations enabled use
#    DEBUGFLAG= -O3 -g -debug extended
    DEBUGFLAG=-g -O0 -warn -traceback -fpe0 -fp-stack-check -gen-interfaces -check all -check noarg_temp_created
    STATICFLAG=
    MATHTRAP=
    EXE=torus.complexity
    STATICLIBS= -Bstatic
    CFITSIOLIBS = -L${CFITSTIOLIBDIR} -lcfitsio
#    CFITSIOLIBS = -L/home/dc-harr1/cfitsio -lcfitsio
    ifeq ($(hdf), yes)
       F95COMP+=-I/cm/shared/apps/hdf5-par/intel/1.6.10/lib
       STARLIBS=-L/cm/shared/apps/hdf5-par/intel/1.6.10/lib
    endif
    ZLIBLIBS = -lz
    KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), davesmac)
#    F95COMP  = g95
#    DEBUGFLAG= -g -fbounds-check -ftrace=full -fimplicit-none -Wall -Wno=140,141,112
    F95COMP  = ifort
#   DEBUGFLAG = -warn -g -no-ipo  # 
    DEBUGFLAG=-O0 -warn -traceback -debug full -check pointer -CB -ftz -save-temps -shared-intel #-assume nounderscore
#    DEBUGFLAG=-g -save-temps -assume nounderscore
#    FASTFLAG = -fast -fp-model fast=2 -no-prec-div -no-prec-sqrt -mdynamic-no-pic -ftz -auto
#    FASTFLAG = -mtune=core2 -xT -O3 -mdynamic-no-pic -no-prec-div -fp-model fast=2 -no-prec-sqrt -shared-intel -prof-genx#-g -save-temps -traceback -mmacosx-version-min=10.4 -intel-static#-ipo #-g
#    FASTFLAG = -mtune=core2 -xT -O3 -mdynamic-no-pic -no-prec-div -fp-model fast=2 -no-prec-sqrt -shared-intel -prof-use #-g -save-temps -traceback -mmacosx-version-min=10.4 -intel-static#-ipo #-g
    FASTFLAG = -mtune=core2 -xT -O3 -mdynamic-no-pic -no-prec-div -fp-model fast=2 -no-prec-sqrt -shared-intel -ip -pad#-g -save-temps -traceback -mmacosx-version-min=10.4 -intel-static#-ipo #-g-ip
#    CFITSIOLIBS = -L/Applications/scisoft/i386/Packages/cfitsio/lib/ -lcfitsio
    CFITSIOLIBS = -L/Users/drundle/bin/cfitsio/lib/ -lcfitsio
    ZLIBLIBS = -lz	
    MKLIBS = -I/Library/Frameworks/Intel_MKL.framework/Versions/10.0.012/Headers -L/Library/Frameworks/Intel_MKL.framework/Versions/10.0.012/Libraries -lmkl_lapack95 -lmkl_intel -lmkl_sequential -lmkl_core -lguide -lpthread #pthread for macs 
#    MKLIBS = -I/Library/Frameworks/Intel_MKL.framework/Versions/10.0.012/Headers -L/Library/Frameworks/Intel_MKL.framework/Versions/10.0.012/Libraries -lmkl_lapack -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lguide -lpthread #pthread for macs 
    EXE = torus.intelmac
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), zen)
    F95COMP= mpif90

    FASTFLAG= -O3 # -ipo
#    FASTFLAG= -O3 -xS #-fp-model fast=2 -no-prec-div -no-prec-sqrt -auto
#    FASTFLAG= -fp-model fast=2 -no-prec-div -no-prec-sqrt -O3 -xT #-ipo 

    DEBUGFLAG=-g -O0 -warn -traceback -fpe0 -fp-stack-check -gen-interfaces -check all -check noarg_temp_created
# Uncomment this line to include MPI correctness checking 
# It is not enabled by default as it can be slow
#    DEBUGFLAG += -check_mpi

## This section was required for setting debug flags with Intel MPI version 3.x 
## There was a conflict between -check in ifort and mpif90. This is fixed in the 4.x versions. 
##
## # Include run-time error checking using ifort. Use environment variable to avoid clash with 
## # ITC -check option in version 3.1 of mpiifort/mpif90 scripts. 
##    ifeq ($(debug), yes)
##      export F90FLAGS=-check all -check noarg_temp_created
##    endif
## # Uncomment this line to include MPI correctness checking using ITC libraries
## # when using version 3.1 of the Intel MPI libraries
## #  DEBUGFLAG += -check

    STATICFLAG=
    MATHTRAP=
    EXE=torus.zen
    STATICLIBS= -Bstatic
    CFITSIOLIBS = -L/home/tjharrie/cfitsio -lcfitsio
    ifeq ($(hdf), yes)
       F95COMP+=-I/sw/hdf5_intel/include
       STARLIBS=-L/sw/hdf5_intel/lib
    endif
    ZLIBLIBS = -lz	
    MKLIBS = -lmkl_lapack95 -lmkl_def -lmkl_em64t -lmkl_intel -lmkl_sequential -lmkl_core -lguide -lpthread #pthread for macs    
    MPIFLAG = -DMPI
    KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), zensingle)
    F95COMP= ifort
    FASTFLAG= -O3 
# For debugging with optimisations enabled use
#    FASTFLAG= -O3 -g -debug extended
    DEBUGFLAG=-g -O0 -warn -traceback -fpe0 -fp-stack-check -gen-interfaces -check all -check noarg_temp_created
    STATICFLAG=
    MATHTRAP=
    EXE=torus.zensingle
    STATICLIBS= -Bstatic
    CFITSIOLIBS = -L/home/tjharrie/cfitsio -lcfitsio
    ifeq ($(hdf), yes)
       F95COMP+=-I/sw/hdf5_intel/include
       STARLIBS=-L/sw/hdf5_intel/lib
    endif
    ZLIBLIBS = -lz
    MPIFLAG = 
    KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), miracle)
    F95COMP  = f90
    DEBUGFLAG= -g
    FASTFLAG = -O -OPT:Olimit=0 -mp -64
#    FASTFLAG = -OPT:Olimit=0 -O -mp
    EXE = torus.miracle
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), ukaff1a)
    F95COMP  = mpif90
    F95FLAG  = -qsuffix=f=f90
    DEBUGFLAG= -g
#    FASTFLAG = -O -OPT:Olimit=0 -mp -64
#    FASTFLAG = -OPT:Olimit=0 -O -mp
    CFITSIOLIBS = -L/home/th/cfitsio -lcfitsio
    EXE = torus.ukaff1a
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

# xlf90 is the IBM compiler
ifeq ($(SYSTEM), omp)
    F95COMP  = xlf90_r
    F95FLAG  = -qsuffix=f=f90 -q64 -qsmp=omp -qnoipa
    DEBUGFLAG= -g
#    FASTFLAG = -O -OPT:Olimit=0 -mp -64
#    FASTFLAG = -OPT:Olimit=0 -O -mp
    EXE = torus.omp
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), ceres)
    F95COMP  = mpiifort
    FASTFLAG  = -O3 -xhost
    DEBUGFLAG=-g -O0 -traceback -fpe0 
    CFITSIOLIBS = -L/home/tjharrie/cfitsio -lcfitsio
    ZLIBLIBS = -lz
    MPIFLAG = -DMPI
    EXE = torus.ceres
    KNOWN_SYSTEM = yes
endif

#
# Deal with make options
#
ETCFLAG  =

ifeq ($(cfitsio), no)
   CFITSIOLIBS  =
else
   ETCFLAG += -DUSECFITSIO
endif

ifeq ($(hdf), yes)
   ETCFLAG += -DUSEHDF
   STARLIBS += -lhdf5_fortran -lhdf5
endif

ifeq ($(zlib), no)
   ZLIBLIBS  =
else
   ETCFLAG += -DUSEZLIB
endif

ifeq ($(memcheck), yes)
   ETCFLAG += -DMEMCHECK
endif

ifeq ($(USEMKL), yes)
   ETCFLAG += -DUSEMKL
else
   MKLIBS = 
endif

ifeq ($(debug), yes)
   FASTFLAG =
   ETCFLAG += $(DEBUGFLAG)
endif

ifeq ($(profile), yes)
   ETCFLAG += -pg
endif

# Coverage flags for gfortran
ifeq ($(coverage), yes)
   ETCFLAG +=  -fprofile-arcs -ftest-coverage
   FASTFLAG = 
endif

# Coverage flags for ifort
ifeq ($(codecov), yes)
   ETCFLAG += -prof-gen=srcpos
   FASTFLAG = 
endif

# Intel Trace Analyser and Collector
ifeq ($(trace), yes)
   TRACEFLAG = -trace
   ETCFLAG  += -tcollect
else
   TRACEFLAG = 
endif

ifeq ($(static), yes)
   ETCFLAG += $(STATICFLAG)
else
   STATICLIBS =
endif

# Set OpenMP flags
ifeq ($(openmp), yes)
  OPENMPFLAG = -openmp
  ifeq ($(SYSTEM), gfortran)
     OPENMPFLAG = -fopenmp 
# -static-libgcc fixes linking problems under Mac OS
     STARLIBS += -static-libgcc
  endif
  ifeq ($(SYSTEM), ompiosx)
     OPENMPFLAG = -fopenmp 
     STARLIBS += -static-libgcc
  endif
  ETCFLAG += $(OPENMPFLAG)
endif

# Add flags required for sphNG compatibility
ifeq ($(withsphng), yes)
  ETCFLAG +=  -mcmodel=medium -i-dynamic
endif

ifeq ($(hydro), yes)
  ETCFLAG += -DHYDRO
endif

ifeq ($(molecular), yes)
  ETCFLAG += -DMOLECULAR
  datacube = yes
endif

ifeq ($(photoion), yes)
  ETCFLAG += -DPHOTOION
  stateq = yes
endif

ifeq ($(xray), yes)
  ETCFLAG += -DXRAY
endif

ifeq ($(sph), yes)
  ETCFLAG += -DSPH
endif

ifeq ($(atomic), yes)
  ETCFLAG += -DCMFATOM
  stateq = yes
  datacube = yes
endif

ifeq ($(stateq), yes)
  ETCFLAG += -DSTATEQ
endif

ifeq ($(datacube), yes)
  ETCFLAG += -DFITSCUBE
endif

ifeq ($(v2), yes)
   TORUSMAINSOURCE= torusMainV2.F90
   INPUTSMODSOURCE= inputs_modV2.F90
else	
   TORUSMAINSOURCE= torusMain.F90
   INPUTSMODSOURCE= inputs_mod.F90
endif

#
# Source files
#

VERSIONSOURCE = torus_version_mod.f90 

F90SOURCE = kind_mod.f90 constants_mod.f90 atom_mod.f90 vector_mod.f90 image_utils_mod.f90 \
            sed_mod.f90 utils_mod.f90  phfit2.f90  bhmie_mod.f90 mieDistCrossSection_mod.f90 \
            gaussian_mod.f90 math_mod.f90 dimensionality_mod.f90 \
            phasematrix_mod.f90 photon_mod.f90 blob_mod.f90 \
            distortion_mod.f90  \
            disc_hydro_mod.f90 gas_opacity_mod.f90 \
            integratePath.f90 fillGridTio.f90 fillGridRayleigh.f90 \
            fillGridThomson.f90 amr_utils_mod.f90  \
            jets_mod.f90 h21cm_mod.f90 \
            source_mod.f90 spectrum_mod.f90 density_mod.f90 \
            sph_data_class.f90 cluster_class.f90 vh1_mod.f90 \
            timing.f90  isochrone_class.f90 \
            filter_set_class.f90 ostar_mod.f90  \
            surface_mod.f90  math_mod2.f90 disc_class.f90 \
            hyd_col_coeff.f90 clump_mod.f90 \
            discwind_class.f90 zeus_data_class.f90 luc_cir3d_class.f90 \
            cmfgen_class.f90 modelatom_mod.f90  magField.f90 \
            romanova_class.f90 messages_mod.f90 starburst_mod.f90 \
            magnetic_mod.f90 units_mod.f90 interferometry_mod.f90 \
            qShep3D_mod.f90 qShep2D_mod.f90 $(VERSIONSOURCE)
            


CAPF90SOURCE =	 amr_mod.F90  octal_mod.F90 mpi_global_mod.F90 random_mod.F90 unix_mod.F90 photoion_utils_mod.F90 photoionAMR_mod.F90 photoion_mod.F90 analytical_velocity_mod.f90 \
	         mpi_amr_mod.F90 hydrodynamics_mod.F90 cmf_mod.F90  nbody_mod.F90 ion_mod.F90 gridtype_mod.F90 gridFromFitsFile.F90 \
	         gridFromFlash.F90 fits_utils_mod.F90 diffusion_mod.F90 formal_solutions.F90 lucy_mod.F90 molecular_mod.F90 \
	         parallel_mod.F90 stateq_mod.F90 image_mod.F90 datacube_mod.F90 vtk_mod.F90 grid_mod.F90 TTauri_mod.F90 angularImage_utils.F90 \
                 phaseloop_mod.F90 gridio_mod.F90 angularImage_mod.F90  timedep_mod.F90 zlib_mod.F90 dust_mod.F90 mieDistPhaseMatrix.F90  \
	         setupamr_mod.F90 physics_mod.F90 outputs_mod.F90 memory_mod.F90 viscosity_mod.F90 xray_mod.F90 $(TORUSMAINSOURCE)

CAPF90SOURCELIB = mpi_global_mod.F90 random_mod.F90 unix_mod.F90  memory_mod.F90 setupamr_mod.F90 physics_mod.F90 outputs_mod.F90 \
	hydrodynamics_mod.F90 angularImage_mod.F90 \
	photoionAMR_mod.F90 photoion_mod.F90 mpi_amr_mod.F90 cmf_mod.F90 \
	          diffusion_mod.F90 formal_solutions.F90 lucy_mod.F90 molecular_mod.F90 \
	          parallel_mod.F90 stateq_mod.F90 torusMod.F90 image_mod.F90 datacube_mod.F90 vtk_mod.F90 \
                  phaseloop_mod.F90 gridio_mod.F90 xray_mod.F90

# Object files
#
INPUTSMODOBJ = $(INPUTSMODSOURCE:.f90=.o)
F90OBJ = $(F90SOURCE:.f90=.o)
CAPF90OBJ =  $(CAPF90SOURCE:.F90=.o)
CAPF90OBJLIB =  $(CAPF90SOURCELIB:.F90=.o)
ALLOBJ = $(CAPF90OBJ) $(F90OBJ)  $(INPUTSMODOBJ)
LIBOBJ = $(CAPF90OBJLIB) $(F90OBJ)  $(INPUTSMODOBJ)

#
# Compilers and compiler flags
#

F90 = $(F95COMP) $(F95FLAG) $(FASTFLAG) $(ETCFLAG) $(MATHTRAP) $(MPI_INC) $(MPIFLAG) $(XMLINCLUDE)

#
#


all: depends $(ALLOBJ)
	svnversion
    ifeq ($(KNOWN_SYSTEM), yes)
	    @echo ""
	    @echo "Linking object files and libraries..."
	    $(F90) $(TRACEFLAG) -o $(EXE) \
	    $(INPUTSMODOBJ) $(F90OBJ) $(CAPF90OBJ) \
	    $(CFITSIOLIBS)  $(STARLIBS) $(STATICLIBS) $(MKLIBS) $(XMLLIBS) $(ZLIBLIBS)
    else
	    @echo ""
	    @echo "make: WARNING: value of SYSTEM = $(SYSTEM) not recognised..."
	    @echo "               ...setting SOURCE_VARIANT to unknown"
    endif

lib: depends $(LIBOBJ)
    ifeq ($(KNOWN_SYSTEM), yes)
	    @echo ""
	    @echo "Linking object files making a library..."
	    $(AR) $(ARFLAGS) libtorus.a $(F90OBJ) $(CAPF90OBJLIB) $(INPUTSMODOBJ)
	    @echo "...done"
    else
	    @echo ""
	    @echo "make: WARNING: value of SYSTEM = $(SYSTEM) not recognised..."
	    @echo "               ...setting SOURCE_VARIANT to unknown"
    endif

preprocess:
	@echo ""
	@echo "====== Creating TORUS for $(SYSTEM) system ======"
	@echo ""
	@echo "BUILD OPTIONS: "
	@echo "OpenMP = $(openmp)"
	@echo "Debug  = $(debug)"


depends: preprocess
	@echo ""
	@echo "Making dependency file..."
	@touch -c make.depends
	@./makedepend
	@echo "...done"
	@echo ""

clean:
	@/bin/rm -f *.o *~ *.mod *.il *.pc *.d *.rnreg *.life *.cfg *.cse *.combine

distclean: clean
	@/bin/rm -f make.depends tags

ctags:
	@ctags --langmap='fortran:+(*.F90)' \
	       *.f90 *.f *.F90

svnver:
    ifeq ($(getsvnver), yes)
	    @echo "calling createsvnversion"
	    @./createsvnversion.csh 
	    @echo "...done"
    endif

torus_version_mod.o: svnver $(VERSIONSOURCE)

link: depends $(ALLOBJ)
    ifeq ($(KNOWN_SYSTEM), yes)
	    @echo ""
	    @echo "Linking object files and libraries..."
	    $(F90) $(TRACEFLAG) -o $(EXE) \
	    $(F90OBJ) $(CAPF90OBJ) \
	    $(CFITSIOLIBS) $(STARLIBS) $(STATICLIBS)
	    @echo "...done"
    else
	    @echo ""
	    @echo "make: WARNING: value of SYSTEM = $(SYSTEM) not recognised..."
	    @echo "               ...setting SOURCE_VARIANT to unknown"
    endif
#
# General rules to make object files
#
.SUFFIXES:

%.o: %.f90
	$(F90) -c $<

%.o: %.F90
	$(F90) -c $<


#
# Auto-generated dependency file
#   - 'make depends' must be run before 'make' as the dependencies included
#     with the directive below are always one step behind the current version
#
-include make.depends

# DO NOT DELETE
