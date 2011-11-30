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
#   10. make zlib=yes:   : Use z-lib compression libary to produce vtu files
#-------------------------------------------------------------------------

USEMKL  = no
cfitsio = yes
memcheck = yes
zlib = yes
debug   = no
profile = no
static  = no
trace   = no
openmp  = no 
v2      = yes
coverage = no
withsphng = no
hydro = yes
molecular = yes
photoion = yes

PATH2G2C     = /usr/lib/gcc-lib/i386-redhat-linux/3.4.4
GUANINE_HOME = /Network/Guanine

KNOWN_SYSTEM = no
MPIFLAG = 
#
# Set system-dependent variables
#

# NAG Fortran compiler
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
#    STARLIBS = -limf                                                                   
    EXE = torus.intelmac
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), g95)
    F95COMP  = g95
    DEBUGFLAG= -g3 -fbounds-check -ftrace=full -fimplicit-none -freal=nan -Wall -Wunused-internal-procs -Wunused-parameter -Wno=140,141,112 -Werror
#    FASTFLAG = -O3 -mcpu=G5
    FASTFLAG = -O3
#    FASTFLAG = -O -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math	
    CFITSIOLIBS = -lcfitsio
    ZLIBLIBS = -lz	
#    F2CLIB   = -L$(PATH2G2C) -lg2c
    EXE = torus.g95
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), sun)
    F95COMP  = sunf95
    DEBUGFLAG= -g
    FASTFLAG = -O3
#    FASTFLAG = -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    F2CLIB   = -L$(PATH2G2C) -lg2c
    EXE = torus.sun
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), osx)
    F95COMP  = g95
    DEBUGFLAG= -g
    FASTFLAG =
#    FASTFLAG = -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    F2CLIB   = -L/usr/local/lib -lg2c
    EXE = torus.osx
    KNOWN_SYSTEM = yes
endif

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
#    F2CLIB   = -L/usr/local/lib -lg2c
    EXE = torus.gfortran
    KNOWN_SYSTEM = yes
    CFITSIOLIBS = -L/Users/${USER}/cfitsio/lib -lcfitsio
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

ifeq ($(SYSTEM), gfort)
    F95COMP  = gfc
    DEBUGFLAG= -g 
#    FASTFLAG = -O3
    FASTFLAG = -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -static
    F2CLIB   = -L$(PATH2G2C) -lg2c
    EXE = torus.g95
    KNOWN_SYSTEM = yes
endif

ifeq  ($(SYSTEM), linux)
    F95COMP  = f95
    #F95FLAG  = -maxcontin=2500
    DEBUGFLAG= -g -C=all
    #DEBUGFLAG= -g
    STATICFLAG = -Bstatic
    FASTFLAG = -w 
    MATHTRAP = -ieee=stop
    F2CLIB   = -lg2c
    STATICLIBS = -lpthread
    EXE = torus.linux
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), alpha)
    F95COMP  = f90
    DEBUGFLAG= -g
    FASTFLAG = -fast
    MATHTRAP = -fpe3
    EXE = torus
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), lamosx)
    F95COMP  = ${GUANINE_HOME}/csr201/local-ppc/bin/mpif77 
    DEBUGFLAG= -g -ftrace=full
    FASTFLAG = -O0
    CFITSIOLIBS = -L${GUANINE_HOME}/tjharrie/cfitsio -lcfitsio
    EXE = torus.lamosx
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), lamosx-intel)
    F95COMP  = ${GUANINE_HOME}/csr201/local-i386/bin/mpif77
    F95FLAG  = -fno-pic -read_only_relocs suppress
    DEBUGFLAG= -g -ftrace=full
    FASTFLAG = -O0
    CFITSIOLIBS = -L${GUANINE_HOME}/tjharrie/cfitsio.i386 -lcfitsio
    EXE = torus.lamosx-intel
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
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

ifeq ($(SYSTEM), ompiosx-intel)
    F95COMP  = ${GUANINE_HOME}/csr201/local-i386/opt/openmpi/bin/mpif77
    F95FLAG  = -fno-pic -read_only_relocs suppress -w
    DEBUGFLAG= -g -ftrace=full
    FASTFLAG = -O0
    CFITSIOLIBS = -L${GUANINE_HOME}/tjharrie/cfitsio.i386 -lcfitsio
    STARLIBS =  -L/usr/lib -lSystemStubs 
    EXE = torus.ompiosx-intel
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), mpilam)
    TOP  = /h/rk/mpi_nag
    TOP2 = /h/rk/mpi_g77
    F95COMP  = ${TOP}/bin/mpif95
    F95FLAG  = -mismatch
    DEBUGFLAG= -g -C=all -dcfuns 
    STATICFLAG = -Bstatic
    FASTFLAG = -dcfuns -w 
#    MATHTRAP = -ieee=stop
    MPI_INC  = -I${TOP}/include
    F2CLIB   = ${MPI_INC} -L${TOP}/lib -lmpi -llamf77mpi -L$(PATH2G2C) -lg2c -lm
    EXE = torus.mpilam
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), mpich)
    TOP  = /h/th/mpich
    TOP2 = /h/th/mpich_g77
    F95COMP  = ${TOP}/bin/mpif90
    F95FLAG  = -mismatch
    DEBUGFLAG= -g -dcfuns
    FASTFLAG = -O3 -dcfuns -w
#    MATHTRAP = -ieee=stop
    MPI_INC  = -I${TOP}/include
    F2CLIB   = ${MPI_INC} -L${TOP}/lib -lmpichf90 -lmpich -L$(PATH2G2C) -lg2c -lm
    EXE = torus.mpich
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), linuxintel)
    F95COMP  = ifort 
    DEBUGFLAG= -g
#    FASTFLAG = -fpp -auto -I/usr/include/nptl -L/usr/lib/nptl
#    FASTFLAG = -g -openmp -fpp -auto -I/usr/include/nptl -L/usr/lib/nptl
#    FASTFLAG = -I/usr/include/nptl -L/usr/lib/nptl -zero 
    CFITSIOLIBS = -lcfitsio
    F2CLIB   = -L/usr/lib -lg2c
    EXE = torus.linuxintel
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), zen)
    F95COMP= mpif90

#    FASTFLAG= -O3 -xS #-fp-model fast=2 -no-prec-div -no-prec-sqrt -auto

    DEBUGFLAG=-g -O0 -warn -traceback -fpe0 -fp-stack-check -gen-interfaces
#    FASTFLAG= -fp-model fast=2 -no-prec-div -no-prec-sqrt -O3 -xT #-ipo 
    FASTFLAG= -O3 # -ipo
#    FASTFLAG= -O3 -xS
#    DEBUGFLAG=-g -O0 -warn -warn nounused -traceback -fpe0 -fp-stack-check -gen-interfaces

# Include run-time error checking using ifort. Use environment variable to avoid clash with 
# ITC -check option in version 3.1 of mpiifort/mpif90 scripts. 
    ifeq ($(debug), yes)
      export F90FLAGS=-check all -check noarg_temp_created
    endif
# Uncomment this line to include MPI correctness checking using ITC libraries
# when using version 3.1 of the Intel MPI libraries
#  DEBUGFLAG += -check 
    STATICFLAG=
    MATHTRAP=
    EXE=torus.zen
    F2CLIB=
#    -L$(PATH2G2C) -lg2c
    STATICLIBS= -Bstatic
    STARLIBS= #-DMPI -Xlinker -rpath -Xlinker $libdir -Xlinker -rpath -Xlinker /opt/intel/mpi-rt/3.0 -lmpi -lmpiif -lmpigi -lrt -lpthread -ldl
    CFITSIOLIBS = -L/home/tjharrie/cfitsio -lcfitsio
    ZLIBLIBS = -lz	
    MKLIBS = -lmkl_lapack95 -lmkl_def -lmkl_em64t -lmkl_intel -lmkl_sequential -lmkl_core -lguide -lpthread #pthread for macs    
    MPIFLAG = -DMPI
    KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), zensingle)
    F95COMP= ifort
    FASTFLAG= -O3 -xT
    DEBUGFLAG=-g -O0 -warn -traceback -fpe0 -fp-stack-check -gen-interfaces -check all -check noarg_temp_created
    STATICFLAG=
    MATHTRAP=
    EXE=torus.zensingle
    F2CLIB=
    STATICLIBS= -Bstatic
    STARLIBS= 
    CFITSIOLIBS = -L/home/tjharrie/cfitsio -lcfitsio
    ZLIBLIBS = -lz
    MPIFLAG = 
    KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), mpiintel)
    TOP  = $(HOME)/mpi_intel
    TOP2 = $(HOME)/mpi_g77
    F95COMP  = ${TOP}/bin/mpif77
    F95FLAG  = -mp
    DEBUGFLAG= -g -O0 -CA -CV -CS -CB -mp
#    FASTFLAG = -O2 -I/usr/include/nptl -L/usr/lib/nptl -mp
    FASTFLAG = -O -mp
    MPI_INC  = -I${TOP}/include
    F2CLIB   = ${MPI_INC} -L${TOP}/lib -lmpi -llamf77mpi -L$(PATH2G2C) -lg2c -lm
    EXE = torus.mpiintel
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), openmpi)
    TOP = /h/chris/local/opt/openmpi
    F95COMP  = ${TOP}/bin/mpif77
    F95FLAG  = -w
    DEBUGFLAG= -g -ftrace=full
#    FASTFLAG = -O3
    MPI_INC  = -I${TOP}/include
    CFITSIOLIBS = -lcfitsio
    F2CLIB   = ${MPI_INC} -L$(PATH2G2C) -lg2c
    EXE = torus.openmpi
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), openmpi-nocona)
    TOP = /home/chris/local
    F95COMP  = ${TOP}/bin/mpif77
    F95FLAG  = -w
    DEBUGFLAG= -g -ftrace=full
#    FASTFLAG = -O0
#    FASTFLAG = -O3
    FASTFLAG = -march=pentium4 -O0
#    FASTFLAG = -march=pentium4 -O3
#    FASTFLAG = -march=nocona -O1
    MPI_INC  = -I${TOP}/include
    CFITSIOLIBS = -lcfitsio
    F2CLIB   = ${MPI_INC} -L$(PATH2G2C) -lg2c
    EXE = torus.openmpi-nocona
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
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

ifeq ($(SYSTEM), mpig95)
    TOP = /h/chris/local
    F95COMP  = ${TOP}/bin/mpif77
    F95FLAG  = -w
    DEBUGFLAG= -g -ftrace=full
#    FASTFLAG = -O3
    MPI_INC  = -I${TOP}/include
    CFITSIOLIBS = -lcfitsio
    F2CLIB   = ${MPI_INC} -L${TOP}/lib -lmpi -llamf77mpi -L$(PATH2G2C) -lg2c -lm
    EXE = torus.mpig95
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), darwin)
    F95COMP  = g95
    DEBUGFLAG= -g -fbounds-check -Wall
    FASTFLAG = -O3
    ZLIBLIBS = -lz	
    CFITSIOLIBS = -L${GUANINE_HOME}/tjharrie/cfitsio -lcfitsio
    STARLIBS =  -L/usr/lib -lSystemStubs
    EXE = torus.darwin
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), mpidarwin)
    TOP = ${HOME}/mpi_g95
    F95COMP  = ${TOP}/bin/mpif77 
    DEBUGFLAG= -g  
    FASTFLAG = -O3
    STARLIBS =  -L/usr/lib -lSystemStubs 
    EXE = latticeQED
    MPIFLAG = -DMPI
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

ifeq ($(coverage), yes)
   ETCFLAG +=  -fprofile-arcs -ftest-coverage
   FASTFLAG = 
endif

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
     FASTFLAG = 
  endif
  ifeq ($(SYSTEM), ompiosx)
     OPENMPFLAG = -fopenmp 
     STARLIBS += -static-libgcc
     FASTFLAG = 
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
endif

ifeq ($(photoion), yes)
  ETCFLAG += -DPHOTOION
endif

ifeq ($(v2), yes)
   TORUSMAINSOURCE= torusMainV2.F90
   INPUTSMODSOURCE= inputs_modV2.f90
else	
   TORUSMAINSOURCE= torusMain.F90
   INPUTSMODSOURCE= inputs_mod.f90
endif

#
# Source files
#

VERSIONSOURCE = torus_version_mod.f90 

F90SOURCE = kind_mod.f90 constants_mod.f90 atom_mod.f90 vector_mod.f90 image_utils_mod.f90 \
            sed_mod.f90 utils_mod.f90  phfit2.f90  bhmie_mod.f90 mieDistCrossSection_mod.f90 \
            gaussian_mod.f90 math_mod.f90 dimensionality_mod.f90 \
            phasematrix_mod.f90 photon_mod.f90 blob_mod.f90 \
            distortion_mod.f90  TTauri_mod.f90 \
            disc_hydro_mod.f90 gas_opacity_mod.f90 \
            integratePath.f90 readIntPro.f90 \
            readTioCrossSection.f90  \
            fillGridTio.f90 fillGridRayleigh.f90 \
            fillGridThomson.f90 amr_utils_mod.f90  \
            jets_mod.f90 h21cm_mod.f90 \
            source_mod.f90 spectrum_mod.f90 \
            wr104_mod.f90 density_mod.f90 \
            sph_data_class.f90 cluster_class.f90 vh1_mod.f90 \
            timing.f90  isochrone_class.f90 \
            filter_set_class.f90 ostar_mod.f90  \
            surface_mod.f90  math_mod2.f90 disc_class.f90 \
            hyd_col_coeff.f90 flow_speed_variables.f90 clump_mod.f90 \
            discwind_class.f90 zeus_data_class.f90 luc_cir3d_class.f90 \
            cmfgen_class.f90  \
            modelatom_mod.f90  magField.f90 \
            romanova_class.f90 messages_mod.f90 starburst_mod.f90 \
            nrtype.f90 nrutil.f90 magnetic_mod.f90 units_mod.f90 interferometry_mod.f90 $(VERSIONSOURCE)
            


CAPF90SOURCE =	 octal_mod.F90 mpi_global_mod.F90 random_mod.F90 unix_mod.F90 photoion_utils_mod.F90 photoionAMR_mod.F90 photoion_mod.F90 \
	         mpi_amr_mod.F90 hydrodynamics_mod.F90 cmf_mod.F90  nbody_mod.F90 ion_mod.F90 gridtype_mod.F90 amr_mod.F90 \
	         fits_utils_mod.F90 diffusion_mod.F90 formal_solutions.F90 lucy_mod.F90 molecular_mod.F90 \
	         parallel_mod.F90 stateq_mod.F90 image_mod.F90 datacube_mod.F90 vtk_mod.F90 grid_mod.F90 \
                 phaseloop_mod.F90 gridio_mod.F90 angularImage_mod.F90  timedep_mod.F90 zlib_mod.F90 dust_mod.F90 mieDistPhaseMatrix.F90  \
	         setupamr_mod.F90 physics_mod.F90 outputs_mod.F90 memory_mod.F90 $(TORUSMAINSOURCE)

CAPF90SOURCELIB = mpi_global_mod.F90 random_mod.F90 unix_mod.F90  memory_mod.F90 setupamr_mod.F90 physics_mod.F90 outputs_mod.F90 \
	hydrodynamics_mod.F90 angularImage_mod.F90 \
	photoionAMR_mod.F90 photoion_mod.F90 mpi_amr_mod.F90 cmf_mod.F90 \
	          diffusion_mod.F90 formal_solutions.F90 lucy_mod.F90 molecular_mod.F90 \
	          parallel_mod.F90 stateq_mod.F90 torusMod.F90 image_mod.F90 datacube_mod.F90 vtk_mod.F90 \
                  phaseloop_mod.F90 gridio_mod.F90 

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
    ifeq ($(KNOWN_SYSTEM), yes)
	    @echo ""
	    @echo "Linking object files and libraries..."
	    $(F90) $(TRACEFLAG) -o $(EXE) \
	    $(INPUTSMODOBJ) $(F90OBJ) $(CAPF90OBJ) \
	    $(CFITSIOLIBS) $(STARLIBS) $(F2CLIB) $(STATICLIBS) $(MKLIBS) $(XMLLIBS) $(ZLIBLIBS)
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
	    @echo "calling createsvnversion"
	    @./createsvnversion.csh 
	    @echo "...done"

torus_version_mod.o: svnver $(VERSIONSOURCE)

link: depends $(ALLOBJ)
    ifeq ($(KNOWN_SYSTEM), yes)
	    @echo ""
	    @echo "Linking object files and libraries..."
	    $(F90) $(TRACEFLAG) -o $(EXE) \
	    $(F90OBJ) $(CAPF90OBJ) \
	    $(CFITSIOLIBS) $(STARLIBS) $(F2CLIB) $(STATICLIBS)
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
