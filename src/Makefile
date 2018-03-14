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
#   12. make sundials=yes: Use sundials chemistry ODE tool
#   13. make f2008=yes   : Use Fortran 2008 features
#-------------------------------------------------------------------------

USEMKL  = no
cfitsio = yes
minusg = no
netcdf = no
sundials = no
hdf = no
f2008=no
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
sph2grid_parallel = no
pdr = no
hydro = yes
molecular = yes
photoion = yes
sph = yes
atomic = yes
xray = yes
chemistry = no
# These are set later on, don't set in the make command
stateq = no
datacube = no
# Run script to get subversion version number
getgitver = yes

KNOWN_SYSTEM = no
MPIFLAG = 

# Assume the default compiler is a recent Intel compiler. OpenMP flags for other compilers
# need to be set in the relevant SYSTEM section of this file
OPENMPFLAG = -qopenmp

# When linking with SPH-NG ifort requires extra flags.
SPHNG_FLAGS =

#
# Set system-dependent variables
#


########################
# NAG Fortran compiler #
########################

ifeq ($(SYSTEM), nagfor)
    F95COMP  = nagfor
    F95FLAG  = -colour -kind=byte -fpp -dcfuns -f2008 -DUSEIEEEISNAN -DUSEUNIXENV
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
    DEBUGFLAG= -ggdb -fbacktrace -fcheck=all -Wunderflow -Wall -Wno-surprising -Wno-conversion -pedantic-errors -std=f2008 -fall-intrinsics -ffpe-trap=invalid,zero,overflow #,underflow,denormal
#    DEBUGFLAG= -ggdb -fbacktrace -fcheck=all -Wunderflow -Wall -Werror -Wno-surprising -pedantic-errors -std=f2008 -fall-intrinsics -ffpe-trap=invalid,zero,overflow #,underflow,denormal
#    DEBUGFLAG= -ggdb -fbacktrace -Wunderflow -Wall -Werror -pedantic-errors -std=f2008 -fall-intrinsics -ffpe-trap=invalid,zero,overflow #,underflow,denormal
    FASTFLAG = -O3
    OPENMPFLAG = -fopenmp
# -static-libgcc fixes linking problems under Mac OS
# -lcurl is required for newer cfitsio libraries to prevent failure at link stage
    STARLIBS += -static-libgcc -static-libgfortran -lcurl
#    FASTFLAG = -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    EXE = torus.gfortran
    KNOWN_SYSTEM = yes
    CFITSIOLIBS = -L${HOME}/cfitsio/lib -lcfitsio
    KROMELIBS = -L${HOME}/krome/build -lkrome
    ifeq ($(chemistry), yes)
       F95COMP+=-I${HOME}/krome/build
    endif

endif

ifeq ($(SYSTEM), thaw)
#    ifeq ($(mpi), yes)
        F95COMP = mpif90
#        F95COMP+=-I/sw/hdf5_intel/include
        MPIFLAG = -DMPI
#    else
#        F95COMP = gfortran
#    endif
     ZLIBLIBS = -lz
#    DEBUGFLAG= -ggdb -fbacktrace -fcheck=all -Wunderflow -Wall -Werror -pedantic-errors -std=f2008 -fall-intrinsics -ffpe-trap=invalid,zero,overflow #,underflow,denormalx
#    DEBUGFLAG= -ggdb -fbacktrace -Wunderflow -Wall -Werror -pedantic-errors -std=f2008 -fall-intrinsics -ffpe-trap=invalid,zero,overflow #,underflow,denormal
    DEBUGFLAG= -ggdb -fbacktrace -Wunderflow -Wall -Werror -pedantic-errors -std=f2008 -fall-intrinsics -ffpe-trap=invalid,zero,overflow #,underflow,denormal


    FASTFLAG = -O3
#    FASTFLAG = -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    EXE = torus.haworth
    KNOWN_SYSTEM = yes
    CFITSIOLIBS = -L/Users/${USER}/cfitsio/lib -lcfitsio
endif

# MPI build with the gfortran compiler. 
ifeq ($(SYSTEM), ompiosx)
    F95COMP  = mpif90 
    DEBUGFLAG= -ggdb -fbacktrace -fcheck=all -Wunderflow -Wall -Werror -Wno-surprising -Wno-error=maybe-uninitialized -pedantic-errors -std=f2008 -fall-intrinsics -ffpe-trap=invalid,zero,overflow #,underflow,denormal
    FASTFLAG = -O3
    OPENMPFLAG = -fopenmp
# -static-libgcc fixes linking problems under Mac OS
    STARLIBS += -static-libgcc
    CFITSIOLIBS = -L${HOME}/cfitsio/lib -lcfitsio
    EXE = torus.ompiosx
    ZLIBLIBS = -lz	
    MPIFLAG = -DMPI
    KROMELIBS = -L${HOME}/krome/build -lkrome
    ifeq ($(chemistry), yes)
       F95COMP+=-I${HOME}/krome/build
    endif
    KNOWN_SYSTEM = yes
endif

# For use with the test suite on post-zen
ifeq ($(SYSTEM), testsuite)
    ifeq ($(mpi), yes)
        F95COMP = mpif90
        MPIFLAG = -DMPI
    else
        F95COMP = gfortran 
    endif
    DEBUGFLAG= -ggdb -fbacktrace -fcheck=all -Wunderflow -Wall -Werror -Wno-surprising -Wno-error=maybe-uninitialized -pedantic-errors -std=f2008 -fall-intrinsics -ffpe-trap=invalid,zero,overflow
    FASTFLAG = -O3
    OPENMPFLAG = -fopenmp
    CFITSIOLIBS = -L/home/torustest/cfitsio/lib -lcfitsio
    EXE = torus.testsuite
    ZLIBLIBS = -lz	
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
    CFITSIOLIBS = -L${HOME}/cfitsio/lib -lcfitsio
    ZLIBLIBS = -lz
    EXE = torus.ifort
    KNOWN_SYSTEM = yes
    SPHNG_FLAGS = -mcmodel=medium
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
    OPENMPFLAG = -fopenmp
    FASTFLAG= -O3 
# For debugging with optimisations enabled use
#    DEBUGFLAG= -O3 -g -debug extended
    DEBUGFLAG=-g -O0 -warn -traceback -fpe0 -fp-stack-check -gen-interfaces -check all -check noarg_temp_created
    STATICFLAG=
    MATHTRAP=
    EXE=torus.complexity
    STATICLIBS= -Bstatic
#    CFITSIOLIBDIR= /home/dc-harr1/cfitsio
#    CFITSIOLIBS = -L${CFITSTIOLIBDIR} -lcfitsio
    CFITSIOLIBS = -L/home/dc-harr1/cfitsio -lcfitsio
    ifeq ($(hdf), yes)
       F95COMP+=-I/cm/shared/apps/hdf5-par/intel/1.6.10/lib
       STARLIBS=-L/cm/shared/apps/hdf5-par/intel/1.6.10/lib
    endif
    KROMELIBS = -L${HOME}/krome/build -lkrome
    ifeq ($(chemistry), yes)
       F95COMP+=-I${HOME}/krome/build
    endif

    SPHNG_FLAGS = -mcmodel=medium -i-dynamic
    ZLIBLIBS = -lz
    KNOWN_SYSTEM=yes
endif

# Isca configuration: to pick up the FITS libraries load the CFITSIO module
# there's no need to add the path to the FITS libraries in this Makefile. 
# The default build uses MPI, for a serial build specify mpi=no. 
ifeq ($(SYSTEM), isca)
    ifeq ($(mpi), no)
       F95COMP= ifort
       MPIFLAG=
    else
       F95COMP= mpiifort
       MPIFLAG= -DMPI
    endif
    FASTFLAG= -O3
# To enable the MPI checking feature in Intel MPI add the -check_mpi flag
# For debugging with optimisations enabled use:
#    DEBUGFLAG= -O3 -g -debug extended
    DEBUGFLAG=-g -O0 -warn -traceback -fpe0 -fp-stack-check -gen-interfaces -check all -check noarg_temp_created
    STATICFLAG=
    EXE=torus.isca
    CFITSIOLIBS = -lcfitsio
    STATICLIBS= -Bstatic
    KROMELIBS = -L${HOME}/krome/build -lkrome
    ifeq ($(chemistry), yes)
       F95COMP+=-I${HOME}/krome/build
    endif
    SPHNG_FLAGS = -mcmodel=medium
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
    F95COMP  = mpif90
    FASTFLAG  = -O3 -xhost
    DEBUGFLAG=-g -O0 -traceback -fpe0 
    CFITSIOLIBS = -L/home/tjharrie/cfitsio -lcfitsio
    ZLIBLIBS = -lz
    MPIFLAG = -DMPI
    EXE = torus.ceres
    SPHNG_FLAGS = -mcmodel=medium -i-dynamic
    KNOWN_SYSTEM = yes
endif


#CFITSIOLIBS = -L/home/thaworth/local/cfitsio/lib -lcfitsio

#
# Deal with make options
#
ETCFLAG  =

ifeq ($(cfitsio), no)
   CFITSIOLIBS  =
else
   ETCFLAG += -DUSECFITSIO
endif

ifeq ($(netcdf), yes)
   ETCFLAG += -DNETCDF
endif

ifeq ($(f2008), yes)
   ETCFLAG += -DF2008
endif

#SUNFLAG = -sundials

ifeq ($(sundials), yes)
#    SUNDIALSLIBRARIES = -L/media/GEN/sundials/v2/lib #-sundials
#    SUNDIALSINC = -I/media/GEN/sundials/v2/include #-sundials
#    SUNDIALSLIBRARIES = -L/media/GEN/sundials/v2/lib #-sundials
#    SUNDIALSINC = -I/media/GEN/sundials/v2/include #-sundials
#    SUNDIALSLIBRARIES = -L/sw/sundials-2.5.0/lib #-sundials                                                                                                                                         
#   SUNDIALSINC = -I/sw/sundials-2.5.0/include #-sundials
    SUNDIALSLIBRARIES = -L/home/th457/sundials/lib #-sundials                                                                                                                                         
    SUNDIALSINC = -I/home/th457/sundials/include #-sundials
    SLIBS = -lsundials_cvode -lsundials_nvecserial -lm
    ETCFLAG += -DUSESUNDIALS
else
#    CSOURCE = 
#    CSOURCELIB =
    SUNDIALSLIBRARIES = 
    SUNDIALSLINC = 
    SLIBS = 
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

ifeq ($(minusg), yes)
   FASTFLAG =
   ETCFLAG += -g
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

# Set OpenMP flag
ifeq ($(openmp), yes)
  ETCFLAG += $(OPENMPFLAG)
endif

# Add flags required for sphNG compatibility
ifeq ($(withsphng), yes)
  ETCFLAG +=  $(SPHNG_FLAGS)
  ETCFLAG +=  -DWITHSPHNG
endif

ifeq ($(sph2grid_parallel), yes)
  ETCFLAG += -DSPH2GRID_PARALLEL
endif

ifeq ($(hydro), yes)
  ETCFLAG += -DHYDRO
endif

ifeq ($(pdr), yes)
  ETCFLAG += -DPDR
  sundials=yes
  CPPFLAG = -cpp
else
  CPPFLAG =
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

ifeq ($(chemistry), no)
   KROMELIBS  =
else
   ETCFLAG += -DCHEMISTRY
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
            distortion_mod.f90   analytical_velocity_mod.f90 \
            disc_hydro_mod.f90 \
            integratePath.f90 fillGridTio.f90 fillGridRayleigh.f90 \
            fillGridThomson.f90 amr_utils_mod.f90  \
            jets_mod.f90 h21cm_mod.f90 \
            source_mod.f90 spectrum_mod.f90 density_mod.f90 \
            cluster_class.f90 vh1_mod.f90 \
            timing.f90  isochrone_class.f90 \
            filter_set_class.f90 ostar_mod.f90  \
            surface_mod.f90  math_mod2.f90 disc_class.f90 \
            hyd_col_coeff.f90 clump_mod.f90 \
            discwind_class.f90 zeus_data_class.f90 luc_cir3d_class.f90 \
            cmfgen_class.f90 modelatom_mod.f90  magField.f90 \
            romanova_class.f90 messages_mod.f90 starburst_mod.f90 \
            magnetic_mod.f90 units_mod.f90 interferometry_mod.f90 \
            qShep3D_mod.f90 qShep2D_mod.f90 turbulence_mod.f90 eos_mod.f90 citations_mod.f90 \
            detector_mod.f90 interface_hierarchy_mod.f90 triangle_mesh_mod.f90 object_media_mod.f90 obj_io_mod.f90 \
            $(VERSIONSOURCE)
            


CAPF90SOURCE =	 amr_mod.F90  octal_mod.F90 mpi_global_mod.F90 random_mod.F90 unix_mod.F90 photoion_utils_mod.F90 photoionAMR_mod.F90 photoion_mod.F90 \
	         mpi_amr_mod.F90 hydrodynamics_mod.F90 cmf_mod.F90  nbody_mod.F90 ion_mod.F90 gridtype_mod.F90 gridFromFitsFile.F90 sph_data_class.F90 \
	         gridFromFlash.F90 fits_utils_mod.F90 diffusion_mod.F90 formal_solutions.F90 lucy_mod.F90 molecular_mod.F90  gas_opacity_mod.F90  \
	         parallel_mod.F90 stateq_mod.F90 image_mod.F90 datacube_mod.F90 vtk_mod.F90 grid_mod.F90 TTauri_mod.F90 angularImage_utils.F90 \
                 phaseloop_mod.F90 gridio_mod.F90 angularImage_mod.F90  timedep_mod.F90 zlib_mod.F90 dust_mod.F90 mieDistPhaseMatrix.F90  \
	         setupamr_mod.F90 physics_mod.F90 outputs_mod.F90 memory_mod.F90 viscosity_mod.F90 xray_mod.F90 loadbalance_mod.F90 \
                 healpix_modules.F90 healpix.F90 nrayshealpix.F90 pdr_mod.F90 pdr_utils_mod.F90 setuppdr_mod.F90 pah_mod.F90 \
		 biophysics_mod.F90 chemistry_mod.F90 gridanalysis_mod.F90 fixing_mod.F90 $(TORUSMAINSOURCE)


CAPF90SOURCELIB = mpi_global_mod.F90 random_mod.F90 unix_mod.F90  memory_mod.F90 setupamr_mod.F90 physics_mod.F90 outputs_mod.F90 \
	hydrodynamics_mod.F90 angularImage_mod.F90 \
	photoionAMR_mod.F90 photoion_mod.F90 mpi_amr_mod.F90 cmf_mod.F90 \
	diffusion_mod.F90 formal_solutions.F90 lucy_mod.F90 molecular_mod.F90 \
	parallel_mod.F90 stateq_mod.F90 torusMod.F90 image_mod.F90 datacube_mod.F90 vtk_mod.F90 \
        phaseloop_mod.F90 gridio_mod.F90 xray_mod.F90 healpix_modules.F90 healpix.F90 nrayshealpix.F90 \
        pdr_mod.F90 pdr_utils_mod.F90 setuppdr_mod.F90 octal_mod.F90 amr_mod.F90 grid_mod.F90 dust_mod.F90 \
        photoion_utils_mod.F90 ion_mod.F90 angularImage_utils.F90 loadbalance_mod.F90 nbody_mod.F90 viscosity_mod.F90 \
        mieDistPhaseMatrix.F90 gas_opacity_mod.F90 TTauri_mod.F90 gridFromFlash.F90 fits_utils_mod.F90 sph_data_class.F90 \
	biophysics_mod.F90 chemistry_mod.F90 gridanalysis_mod.F90 fixing_mod.F90 pah_mod.F90 


ifeq ($(pdr),yes)
CSOURCE = calculate_abundances_onepart.c jacobian_reduced.c odes_reduced.c
CSOURCELIB = calculate_abundances_onepart.c jacobian_reduced.c odes_reduced.c
endif

# Object files
#
#ifeq($(sundials), no)
#    CSOURCE = ##
#    CSOURCELIB =
#endif

INPUTSMODOBJ = $(INPUTSMODSOURCE:.F90=.o)
F90OBJ = $(F90SOURCE:.f90=.o)
CAPF90OBJ =  $(CAPF90SOURCE:.F90=.o)
CAPF90OBJLIB =  $(CAPF90SOURCELIB:.F90=.o)
COBJ = $(CSOURCE:.c=.o)
CLIB = $(CSOURCELIB:.c=.o)
ALLOBJ = $(CAPF90OBJ) $(F90OBJ)  $(INPUTSMODOBJ) $(COBJ)
LIBOBJ = $(CAPF90OBJLIB) $(F90OBJ)  $(INPUTSMODOBJ) $(CLIB)

#
# Compilers and compiler flags
#

F90 = $(F95COMP) $(F95FLAG) $(FASTFLAG) $(ETCFLAG) $(MATHTRAP) $(MPI_INC) $(MPIFLAG) $(XMLINCLUDE) $(NCDF_INC)
GCC = gcc -fopenmp -O0 -Wall $(SUNDIALSINC)
#
#


all: depends $(ALLOBJ)
    ifeq ($(KNOWN_SYSTEM), yes)
	    @echo ""
	    @echo "Linking object files and libraries..."
#	    ifeq($(sundials), yes)
#	    $(GCC) -o  $(COBJ) 
#	    endif
	    $(F90) $(TRACEFLAG) $(SUNDIALSLIBRARIES) -o $(EXE) \
	    $(INPUTSMODOBJ) $(F90OBJ) $(CAPF90OBJ) \
	    $(CFITSIOLIBS)  $(STARLIBS) $(STATICLIBS) $(MKLIBS) $(XMLLIBS) $(ZLIBLIBS) $(COBJ) \
	    $(SLIBS) $(NCDF_LIB) $(KROMELIBS)

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
	@echo "Sundials  = $(sundials)"


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

gitver:
    ifeq ($(getgitver), yes)
	    @echo "calling creategitversion"
	    @./creategitversion 
	    @echo "...done"
    endif

torus_version_mod.o: gitver $(VERSIONSOURCE)

link: depends $(ALLOBJ)
    ifeq ($(KNOWN_SYSTEM), yes)
	    @echo ""
	    @echo "Linking object files and libraries..."
	    $(F90) $(TRACEFLAG) $(SUNDIALSLIBRARIES) -o $(EXE) \
	    $(F90OBJ) $(CAPF90OBJ) \
	    $(CFITSIOLIBS) $(STARLIBS) $(STATICLIBS) $(SLIBS)
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
	$(F90) $(CPPFLAG) -c $<

%.o: %.F90
	$(F90) $(CPPFLAG) -c $<

%.o: %.c
	$(GCC) -c $<

#
# Auto-generated dependency file
#   - 'make depends' must be run before 'make' as the dependencies included
#     with the directive below are always one step behind the current version
#
-include make.depends

# DO NOT DELETE
