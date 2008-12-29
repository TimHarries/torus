#----------------------------------------------------------------
# Usage :: 
#    1. make             : Make torus with optimization flag(s)
#    2. make clean       : Clean up directory
#    3. make cfitsio=no  : Use cfitsio stubs rather than library
#    4. make debug=yes   : Make torus with debug option,
#    5. make profile=yes : Make torus with profiler option
#    6. make static=yes  : Make torus with static linking
#-----------------------------------------------------------------

USEMKL  = no
cfitsio = yes
debug   = no
profile = no
static  = no
trace   = no

PATH2G2C     = /usr/lib/gcc-lib/i386-redhat-linux/3.4.4
GUANINE_HOME = /Network/Guanine

KNOWN_SYSTEM = no
MPIFLAG = 
#
# Set system-dependent variables
#
ifeq ($(SYSTEM), intelmac)
    F95COMP  = g95
    F77COMP  = g95
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g -fbounds-check -ftrace=full -fimplicit-none -Wall -Wno=140,141,112
    STATICFLAG = -static
    FASTFLAG = -O3
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
    CFITSIOLIBS = -L/Users/${USER}/cfitsio/lib -lcfitsio
    EXE = torus.intelmac
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), davesmac)
#    F95COMP  = g95
#    F77COMP  = g95
#    F77FLAG  = -w -fno-second-underscore -g
#    DEBUGFLAG= -g -fbounds-check -ftrace=full -fimplicit-none -Wall -Wno=140,141,112
    F95COMP  = ifort #-I/Library/Frameworks/Intel_MKL.framework/Versions/10.0.012/Headers -L/Library/Frameworks/Intel_MKL.framework/Versions/10.0.012/Libraries
    F77COMP  = ifort
    F77FLAG  = -w -fno-second-underscore -g -O0 -traceback -fpe0 -debug full -check all -warn all -intconstant -f77rtl
#   DEBUGFLAG = -warn -g -no-ipo  # 
    DEBUGFLAG=-g -O0 -warn nounused -warn -traceback -debug full -check pointer -CB -ftz
#    FASTFLAG = -fast -fp-model fast=2 -no-prec-div -no-prec-sqrt -mdynamic-no-pic -ftz -auto
    FASTFLAG = -mtune=core2 -xT -O3 -mdynamic-no-pic -no-prec-div -fp-model fast=2 -no-prec-sqrt -g
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
    CFITSIOLIBS = -L/scisoft/i386/Packages/cfitsio/lib -lcfitsio# -lmkl_intel_thread -lmkl_core -lguide #-lmkl_intel
#    STARLIBS = -L/star/lib `ndf_link`                                                                   
#    F2CLIB   = -L$(PATH2G2C) -lg2c                                                                      
    EXE = torus.intelmac
    KNOWN_SYSTEM = yes
    USEMKL  = no
endif



ifeq ($(SYSTEM), g95)
    F95COMP  = g95
    F77COMP  = g95
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g -fbounds-check -ftrace=full  -fimplicit-none -Wall -Wno=140,141,112 -Werror
#    FASTFLAG = -O3 -mcpu=G5
    FASTFLAG = -O3
#    FASTFLAG = -O -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math	
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    CFITSIOLIBS = -lcfitsio
    STARLIBS = -L/star/lib `ndf_link`
#    F2CLIB   = -L$(PATH2G2C) -lg2c
    EXE = torus.g95
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), sun)
    F95COMP  = sunf95
    F77COMP  = sunf95
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g
    FASTFLAG = -O3
#    FASTFLAG = -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    STARLIBS = -L/star/lib `ndf_link`
    F2CLIB   = -L$(PATH2G2C) -lg2c
    EXE = torus.sun
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), osx)
    F95COMP  = g95
    F77COMP  = g95
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g
    FASTFLAG =
#    FASTFLAG = -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
    F2CLIB   = -L/usr/local/lib -lg2c
    EXE = torus.osx
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), gfortran)
    F95COMP  = gfortran
    F77COMP  = gfortran
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g
    FASTFLAG = -O3
#    FASTFLAG = -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    STARLIBS = -L/star/lib `ndf_link`
    F2CLIB   = -L/usr/local/lib -lg2c
    EXE = torus.gfortran
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), macmpi)
    F95COMP  = g95
    F77COMP  = g95
    F95FLAG  = -fno-second-underscore
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g
    STATICFLAG = -static
#    FASTFLAG = -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    MPI_INC  = -I/Users/th/macmpi
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    EXE = torus.macmpi
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), gfort)
    F95COMP  = gfc
    F77COMP  = gfc
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g 
#    FASTFLAG = -O3
    FASTFLAG = -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -static
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    STARLIBS = -L/star/lib `ndf_link` 
    F2CLIB   = -L$(PATH2G2C) -lg2c
    EXE = torus.g95
    KNOWN_SYSTEM = yes
endif

ifeq  ($(SYSTEM), linux)
    F95COMP  = f95
    F77COMP  = f77
    #F95FLAG  = -maxcontin=2500
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g -C=all
    #DEBUGFLAG= -g
    STATICFLAG = -Bstatic
    FASTFLAG = -w 
    MATHTRAP = -ieee=stop
    UNIXMOD  = unix_mod.f90.linux
    KINDMOD  = kind_mod.f90.linux
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    STARLIBS = -L/star/lib `ndf_link`
    F2CLIB   = -lg2c
    STATICLIBS = -lpthread
    EXE = torus.linux
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), alpha)
    F95COMP  = f90
    F77COMP  = f77
    DEBUGFLAG= -g
    FASTFLAG = -fast
    MATHTRAP = -fpe3
    UNIXMOD  = unix_mod.f90.osf
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    KINDMOD  = kind_mod.f90.osf
    STARLIBS = -L/star/lib `ndf_link`
    EXE = torus
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), lamosx)
    F77COMP  = ${GUANINE_HOME}/csr201/local-ppc/bin/mpif77 
    F95COMP  = ${GUANINE_HOME}/csr201/local-ppc/bin/mpif77 
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g -ftrace=full
    FASTFLAG = -O0
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
    CFITSIOLIBS = -L${GUANINE_HOME}/tjharrie/cfitsio -lcfitsio
    EXE = torus.lamosx
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), lamosx-intel)
    F95COMP  = ${GUANINE_HOME}/csr201/local-i386/bin/mpif77
    F77COMP  = ${GUANINE_HOME}/csr201/local-i386/bin/mpif77
    F95FLAG  = -fno-pic -read_only_relocs suppress
    F77FLAG  = -fno-pic -read_only_relocs suppress -w -fno-second-underscore
    DEBUGFLAG= -g -ftrace=full
    FASTFLAG = -O0
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
    CFITSIOLIBS = -L${GUANINE_HOME}/tjharrie/cfitsio.i386 -lcfitsio
    EXE = torus.lamosx-intel
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), ompiosx)
    F77COMP  = ${GUANINE_HOME}/csr201/local-ppc/opt/openmpi/bin/mpif77 
    F95COMP  = ${GUANINE_HOME}/csr201/local-ppc/opt/openmpi/bin/mpif77 
    F95FLAG  = -w
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g -ftrace=full
    FASTFLAG = -O0
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
    CFITSIOLIBS = -L${GUANINE_HOME}/tjharrie/cfitsio -lcfitsio
    STARLIBS =  -L/usr/lib -lSystemStubs 
    EXE = torus.ompiosx
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), ompi)
    F77COMP  = mpif90
    F95COMP  = mpif90
    F95FLAG  = -Wno=155
    F77FLAG  = -Wno=151
    DEBUGFLAG= -g -ftrace=full -fbounds-check -fimplicit-none -Wall -Wno=140,141,112
    FASTFLAG = -O3 -march=pentium4 -mtune=pentium4 -malign-double
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
    CFITSIOLIBS =  -L/Users/${USER}/cfitsio/lib -lcfitsio
    EXE = torus.ompi
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), ompiosx-intel)
    F95COMP  = ${GUANINE_HOME}/csr201/local-i386/opt/openmpi/bin/mpif77
    F77COMP  = ${GUANINE_HOME}/csr201/local-i386/opt/openmpi/bin/mpif77
    F95FLAG  = -fno-pic -read_only_relocs suppress -w
    F77FLAG  = -fno-pic -read_only_relocs suppress -w -fno-second-underscore
    DEBUGFLAG= -g -ftrace=full
    FASTFLAG = -O0
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
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
    F77COMP  = ${TOP2}/bin/mpif77
    F95FLAG  = -mismatch
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g -C=all -dcfuns 
    STATICFLAG = -Bstatic
    FASTFLAG = -dcfuns -w 
#    MATHTRAP = -ieee=stop
    MPI_INC  = -I${TOP}/include
    UNIXMOD  = unix_mod.f90.linux
    KINDMOD  = kind_mod.f90.linux
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    STARLIBS = -L/star/lib `ndf_link` 
    F2CLIB   = ${MPI_INC} -L${TOP}/lib -lmpi -llamf77mpi -L$(PATH2G2C) -lg2c -lm
    EXE = torus.mpilam
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), mpich)
    TOP  = /h/th/mpich
    TOP2 = /h/th/mpich_g77
    F95COMP  = ${TOP}/bin/mpif90
    F77COMP  = ${TOP2}/bin/mpif77
    F95FLAG  = -mismatch
    DEBUGFLAG= -g -dcfuns
    FASTFLAG = -O3 -dcfuns -w
#    MATHTRAP = -ieee=stop
    MPI_INC  = -I${TOP}/include
    F77FLAG  = -w -fno-second-underscore
    UNIXMOD  = unix_mod.f90.linux
    KINDMOD  = kind_mod.f90.linux
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    STARLIBS = -L/star/lib `ndf_link`
    F2CLIB   = ${MPI_INC} -L${TOP}/lib -lmpichf90 -lmpich -L$(PATH2G2C) -lg2c -lm
    EXE = torus.mpich
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), linuxintel)
    F95COMP  = ifort 
    F77COMP  = ifort
    F77FLAG  = -w
    DEBUGFLAG= -g
#    FASTFLAG = -fpp -auto -I/usr/include/nptl -L/usr/lib/nptl
#    FASTFLAG = -g -openmp -fpp -auto -I/usr/include/nptl -L/usr/lib/nptl
#    FASTFLAG = -I/usr/include/nptl -L/usr/lib/nptl -zero 
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    CFITSIOLIBS = -lcfitsio
    STARLIBS = -L/star/lib `ndf_link` 
    F2CLIB   = -L/usr/lib -lg2c
    EXE = torus.linuxintel
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), zen)
    F95COMP= mpif90
    F77COMP= mpif90
    FASTFLAG= -O3 -xT -mcmodel=medium -i-dynamic
    DEBUGFLAG=-g -O0 -mcmodel=medium -warn -warn nounused -traceback -fpe0 -fp-stack-check -gen-interfaces
# Include run-time error checking using ifort. Use environment variable to avoid clash with 
# ITC -check option in version 3.1 of mpiifort/mpif90 scripts. 
    ifeq ($(debug), yes)
      export F90FLAGS=-check all -check noarg_temp_created
    endif
# Uncomment this line to include MPI correctness checking using ITC libraries
# when using version 3.1 of the Intel MPI libraries
#    DEBUGFLAG += -check 
    STATICFLAG=
    MATHTRAP=
    UNIXMOD=unix_mod.f90.linuxintel
    KINDMOD=kind_mod.f90.linuxintel
    EXE=torus.zen
    F77FLAG=
    IMAGEF77=image_f77.f.miracle
    SPECF77=spec_f77.f.miracle
    F2CLIB=
#    -L$(PATH2G2C) -lg2c
    STATICLIBS= -Bstatic
    STARLIBS= #-DMPI -Xlinker -rpath -Xlinker $libdir -Xlinker -rpath -Xlinker /opt/intel/mpi-rt/3.0 -lmpi -lmpiif -lmpigi -lrt -lpthread -ldl
#    -L/star/lib `ndf_link` 
    CFITSIOLIBS = -L/home/tjharrie/cfitsio -lcfitsio
    MPIFLAG = -DMPI
    KNOWN_SYSTEM=yes
endif

# Build using SGI's Message Passing Toolkit version of MPI. 
 ifeq ($(SYSTEM), zenmpt)
     F95COMP= ifort -lmpi
     F77COMP= ifort -lmpi
     FASTFLAG= -O3 -xT -mcmodel=medium -i-dynamic
     DEBUGFLAG=-g -O0 -mcmodel=medium -warn -warn nounused -traceback -fpe0 -fp-stack-check -gen-interfaces
 # Include run-time error checking using ifort. Use environment variable to avoid clash with 
 # ITC -check option in version 3.1 of mpiifort/mpif90 scripts. 
     ifeq ($(debug), yes)
       export F90FLAGS=-check all -check noarg_temp_created
     endif
     STATICFLAG=
     MATHTRAP=
     UNIXMOD=unix_mod.f90.linuxintel
     KINDMOD=kind_mod.f90.linuxintel
     EXE=torus.mpt
     F77FLAG=
     IMAGEF77=image_f77.f.miracle
     SPECF77=spec_f77.f.miracle
     F2CLIB=
     STATICLIBS= -Bstatic
     STARLIBS= 
     CFITSIOLIBS = -L/home/tjharrie/cfitsio -lcfitsio
     MPIFLAG = -DMPI
     KNOWN_SYSTEM=yes
 endif


ifeq ($(SYSTEM), mpiintel)
    TOP  = $(HOME)/mpi_intel
    TOP2 = $(HOME)/mpi_g77
    F95COMP  = ${TOP}/bin/mpif77
    F77COMP  = ${TOP2}/bin/mpif77
    F95FLAG  = -mp
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g -O0 -CA -CV -CS -CB -mp
#    FASTFLAG = -O2 -I/usr/include/nptl -L/usr/lib/nptl -mp
    FASTFLAG = -O -mp
    MPI_INC  = -I${TOP}/include
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    STARLIBS = -L/star/lib `ndf_link` 
    F2CLIB   = ${MPI_INC} -L${TOP}/lib -lmpi -llamf77mpi -L$(PATH2G2C) -lg2c -lm
    EXE = torus.mpiintel
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), openmpi)
    TOP = /h/chris/local/opt/openmpi
    F95COMP  = ${TOP}/bin/mpif77
    F77COMP  = ${TOP}/bin/mpif77
    F95FLAG  = -w
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g -ftrace=full
#    FASTFLAG = -O3
    MPI_INC  = -I${TOP}/include
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    CFITSIOLIBS = -lcfitsio
    STARLIBS = -L/star/lib `ndf_link` 
    F2CLIB   = ${MPI_INC} -L$(PATH2G2C) -lg2c
    EXE = torus.openmpi
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), openmpi-nocona)
    TOP = /home/chris/local
    F95COMP  = ${TOP}/bin/mpif77
    F77COMP  = ${TOP}/bin/mpif77
    F95FLAG  = -w
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g -ftrace=full
#    FASTFLAG = -O0
#    FASTFLAG = -O3
    FASTFLAG = -march=pentium4 -O0
#    FASTFLAG = -march=pentium4 -O3
#    FASTFLAG = -march=nocona -O1
    MPI_INC  = -I${TOP}/include
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    CFITSIOLIBS = -lcfitsio
    STARLIBS = -L/star/lib `ndf_link` 
    F2CLIB   = ${MPI_INC} -L$(PATH2G2C) -lg2c
    EXE = torus.openmpi-nocona
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), miracle)
    F95COMP  = f90
    F77COMP  = f77
    DEBUGFLAG= -g
    FASTFLAG = -O -OPT:Olimit=0 -mp -64
#    FASTFLAG = -OPT:Olimit=0 -O -mp
    UNIXMOD  = unix_mod.f90.miracle
    KINDMOD  = kind_mod.f90.miracle
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
    EXE = torus.miracle
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), ukaff1a)
    F95COMP  = mpif90
    F77COMP  = mpif77
    F95FLAG  = -qsuffix=f=f90
    DEBUGFLAG= -g
#    FASTFLAG = -O -OPT:Olimit=0 -mp -64
#    FASTFLAG = -OPT:Olimit=0 -O -mp
    UNIXMOD  = unix_mod.f90.miracle
    KINDMOD  = kind_mod.f90.miracle
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
    CFITSIOLIBS = -L/home/th/cfitsio -lcfitsio
    EXE = torus.ukaff1a
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), omp)
    F95COMP  = xlf90_r
    F77COMP  = xlf_r
    F95FLAG  = -qsuffix=f=f90 -q64 -qsmp=omp -qnoipa
    F77FLAG  = -q64 -qsmp=noauto
    DEBUGFLAG= -g
#    FASTFLAG = -O -OPT:Olimit=0 -mp -64
#    FASTFLAG = -OPT:Olimit=0 -O -mp
    UNIXMOD  = unix_mod.f90.miracle
    KINDMOD  = kind_mod.f90.miracle
    SPECF77  = spec_f77.f.miracle
    IMAGEF77 = image_f77.f.miracle
    EXE = torus.omp
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), ceres)
    F95COMP  = xlf95
    F77COMP  = xlf
    F95FLAG  = -bloadmap:outp -qmoddir=/lscratch -I/lscratch -qsuffix=f=f90 -qnolm -O3 -qarch=pwr3 -qmaxmem=-1
    F77FLAG  = -qmoddir=/lscratch -I/lscratch -qsuffix=f=f -qnolm -O3 -qarch=pwr3 -qmaxmem=-1
    DEBUGFLAG  = -g
    UNIXMOD  = unix_mod.f90.aix
    KINDMOD  = kind_mod.f90.aix
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
    EXE = torus.aix
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), mpig95)
    TOP = /h/chris/local
    F95COMP  = ${TOP}/bin/mpif77
    F77COMP  = ${TOP}/bin/mpif77
    F95FLAG  = -w
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g -ftrace=full
#    FASTFLAG = -O3
    MPI_INC  = -I${TOP}/include
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    CFITSIOLIBS = -lcfitsio
    STARLIBS = -L/star/lib `ndf_link` 
    F2CLIB   = ${MPI_INC} -L${TOP}/lib -lmpi -llamf77mpi -L$(PATH2G2C) -lg2c -lm
    EXE = torus.mpig95
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), darwin)
    F95COMP  = g95
    F77COMP  = g95
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g -fbounds-check -Wall
    FASTFLAG = -O3
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
    CFITSIOLIBS = -L${GUANINE_HOME}/tjharrie/cfitsio -lcfitsio
    STARLIBS =  -L/usr/lib -lSystemStubs
    EXE = torus.darwin
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), mpidarwin)
    TOP = ${HOME}/mpi_g95
    F95COMP  = ${TOP}/bin/mpif77 
    F77COMP  = ${TOP}/bin/mpif77 
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g  
    FASTFLAG = -O3
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
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

ifeq ($(USEMKL), yes)
   ETCFLAG += -DUSEMKL
endif

ifeq ($(debug), yes)
   FASTFLAG =
   ETCFLAG += $(DEBUGFLAG)
endif

ifeq ($(profile), yes)
   ETCFLAG += -pg
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

#
# Source files
#
F90SOURCE = kind_mod.f90 constants_mod.f90 atom_mod.f90 vector_mod.f90  \
            unix_mod.f90 inputs_mod.f90 utils_mod.f90 bhmie_mod.f90 mieDistCrossSection_mod.f90 \
            gaussian_mod.f90 grid_mod.f90 math_mod.f90 \
            phasematrix_mod.f90 photon_mod.f90 blob_mod.f90 \
            distortion_mod.f90  TTauri_mod.f90 puls_mod.f90 \
            disc_hydro_mod.f90 gas_opacity_mod.f90 \
            integratePath.f90 readIntPro.f90 mieDistPhaseMatrix.f90 \
            gridtype_mod.f90 getRefractiveIndex.f90 draineCrossSection.f90 \
            readTioCrossSection.f90 rayleighCrossSection.f90 \
            fillGridTio.f90 fillGridRayleigh.f90 fillGridMie.f90 \
            fillGridThomson.f90 cont_read.f90 \
            getMeanMass.f90 input_variables.f90 amr_mod.f90 octal_mod.f90 \
            parameters_mod.f90 jets_mod.f90 \
            dust_mod.f90 source_mod.f90 spectrum_mod.f90 ion_mod.f90 \
            wr104_mod.f90 opacity_lte_mod.f90 density_mod.f90 \
            sph_data_class.f90 cluster_class.f90 \
            linked_list_class.f90 timing.f90  isochrone_class.f90 \
            filter_set_class.f90 ostar_mod.f90  \
            cluster_utils.f90 surface_mod.f90  math_mod2.f90 disc_class.f90 \
            hyd_col_coeff.f90 rotation_variables.f90 \
            flow_speed_variables.f90 clump_mod.f90 \
            discwind_class.f90 zeus_data_class.f90 luc_cir3d_class.f90 \
            jet_class.f90 cmfgen_class.f90  irreg_grid_data_class.f90 \
            modelatom_mod.f90  magField.f90 \
            romanova_class.f90 messages_mod.f90 starburst_mod.f90 \
            mpi_global_mod.f90 nrtype.f90 nrutil.f90 gridio_mod.f90 benchmark_mod.f90 bitstring_mod.f90 

F77SOURCE = phfit2.f spec_f77.f image_f77.f 


CAPF90SOURCE =	 photoionAMR_mod.F90 photoion_mod.F90 mpi_amr_mod.F90 hydrodynamics_mod.F90 cmf_mod.F90 \
	         diffusion_mod.F90 formal_solutions.F90 lucy_mod.F90 molecular_mod.F90 \
	         parallel_mod.F90 stateq_mod.F90 torusMain.F90 image_mod.F90 datacube_mod.F90 vtk_mod.F90 \
                 phaseloop_mod.F90

CAPF90SOURCELIB = photoionAMR_mod.F90 photoion_mod.F90 mpi_amr_mod.F90 cmf_mod.F90 \
	          diffusion_mod.F90 formal_solutions.F90 lucy_mod.F90 molecular_mod.F90 \
	          parallel_mod.F90 stateq_mod.F90 torusMod.F90 image_mod.F90 datacube_mod.F90 vtk_mod.F90 \
                  phaseloop_mod.F90

# Object files
#
F90OBJ = $(F90SOURCE:.f90=.o)
CAPF90OBJ =  $(CAPF90SOURCE:.F90=.o)
CAPF90OBJLIB =  $(CAPF90SOURCELIB:.F90=.o)
F77OBJ = $(F77SOURCE:.f=.o)
ALLOBJ = $(CAPF90OBJ) $(F90OBJ)  $(F77OBJ)
LIBOBJ = $(CAPF90OBJLIB) $(F90OBJ)  $(F77OBJ)

#
# Compilers and compiler flags
#

F90 = $(F95COMP) $(F95FLAG) $(ETCFLAG) $(FASTFLAG) $(MATHTRAP) $(MPI_INC) $(MPIFLAG) $(SPHFLAG)
F77 = $(F77COMP) $(ETCFLAG) $(FASTFLAG) $(MATHTRAP) $(F77FLAG) $(MPI_INC)

#
# Special targets
#   - New .f90.raw files need to be added to both preprocess *and* distclean.
#
all: depends $(ALLOBJ)
    ifeq ($(KNOWN_SYSTEM), yes)
	    @echo ""
	    @echo "Linking object files and libraries..."
	    $(F90) $(TRACEFLAG) -o $(EXE) \
	    $(F77OBJ) $(F90OBJ) $(CAPF90OBJ) \
	    $(CFITSIOLIBS) $(STARLIBS) $(F2CLIB) $(STATICLIBS)
	    @echo "...done"
    else
	    @echo ""
	    @echo "make: WARNING: value of SYSTEM = $(SYSTEM) not recognised..."
	    @echo "               ...setting SOURCE_VARIANT to unknown"
    endif

lib: depends $(LIBOBJ)
    ifeq ($(KNOWN_SYSTEM), yes)
	    @echo ""
	    @echo "Linking object files making a library..."
	    $(AR) $(ARFLAGS) libtorus.a $(F90OBJ) $(CAPF90OBJLIB) $(F77OBJ)
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
	@echo "Creating links for system dependent files..."
	ln -sf $(UNIXMOD) unix_mod.f90
	ln -sf $(KINDMOD) kind_mod.f90
	ln -sf $(IMAGEF77) image_f77.f
	ln -sf $(SPECF77) spec_f77.f
	@echo "...done"

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
	@/bin/rm -f unix_mod.f90 kind_mod.f90 image_f77.f spec_f77.f

ctags:
	@ctags --langmap='fortran:+(*.F90)' \
	       *.f90 *.f *.F90

link: depends $(ALLOBJ)
    ifeq ($(KNOWN_SYSTEM), yes)
	    @echo ""
	    @echo "Linking object files and libraries..."
	    $(F90) $(TRACEFLAG) -o $(EXE) \
	    $(F77OBJ) $(F90OBJ) $(CAPF90OBJ) \
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

%.o: %.f90.raw
	$(F90) -c $*.f90

%.o: %.f90
	$(F90) -c $<

%.o: %.F90
	$(F90) -c $<

%.o: %.f
	$(F77) -c $<

#
# Auto-generated dependency file
#   - 'make depends' must be run before 'make' as the dependencies included
#     with the directive below are always one step behind the current version
#
-include make.depends

# DO NOT DELETE
