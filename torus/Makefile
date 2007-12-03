#----------------------------------------------------------------
# Usage :: 
#    1. make             : Make torus with optimization flag(s)
#    2. make clean       : Clean up directory
#    3. make pgplot=no   : Use pgplot stubs rather than library
#    4. make cfitsio=no  : Use cfitsio stubs rather than library
#    5. make debug=yes   : Make torus with debug option,
#    6. make profile=yes : Make torus with profiler option
#    7. make static=yes  : Make torus with static linking
#-----------------------------------------------------------------

pgplot  = yes
cfitsio = yes
debug   = no
profile = no
static  = no

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
    F77FLAG  = -fno-second-underscore
    DEBUGFLAG= -g -fbounds-check -ftrace=full -fimplicit-none -Wall -Wno=140,141,112
    FASTFLAG = -O3
#    FASTFLAG = -O -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
#   PGLIBS   = -lgfortran -L/usr/local -lpgplot -L/usr/X11R6/lib -lX11 -L/scisoft/i386/lib -lpng12 -laquaterm /System/Library/Frameworks/Cocoa.framework/Cocoa
#   PGLIBS   = -L/usr/local/lib -lpgplot -L/usr/X11R6/lib -lX11 -lpng12 -laquaterm /System/Library/Frameworks/Cocoa.framework/Cocoa
#   PGLIBS   = -L/usr/local/lib /usr/local/lib/libpgplot.a -L/usr/X11R6/lib -lX11 -lpng12 -laquaterm /System/Library/Frameworks/Cocoa.framework/Cocoa
#    CFITSIOLIBS = -L/star/lib/ -lcfitsio
#    PGLIBS   = -L/sw/lib/pgplot -lgfortran -lpgplot -L/usr/X11R6/lib -lX11 -L/sw/lib -lpng -laquaterm /System/Library/Frameworks/Cocoa.framework/Cocoa
#    CFITSIOLIBS = -L/Users/${USER}/cfitsio/lib -lcfitsio
    PGLIBS   = -L/sw/lib/pgplot -lpgplot -L/usr/X11R6/lib -lX11 -L/sw/lib -lpng -laquaterm /System/Library/Frameworks/Cocoa.framework/Cocoa
    CFITSIOLIBS = -L/Users/${USER}/cfitsio/lib -lcfitsio
#    STARLIBS = -L/star/lib `ndf_link`
#    F2CLIB   = -L$(PATH2G2C) -lg2c
    EXE = torus.intelmac
    MK_MPI_SRC = no
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), davesmac)
    F95COMP  = ifort -L/Library/Frameworks/Intel_MKL.framework/Versions/10.0.012/Libraries -I/Library/Frameworks/Intel_MKL.framework/Versions/10.0.012/Headers
    F77COMP  = ifort
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG = -warn -g -no-ipo -debug full# -fbounds-check -ftrace=full                                      
    FASTFLAG = -fast -fp-model fast=2 -no-prec-div -no-prec-sqrt -mdynamic-no-pic -ftz -auto -ipo-jobs2
#    FASTFLAG = -parallel -xT -O3 -ip -fp-model fast=2 -no-prec-div -no-prec-sqrt -mdynamic-no-pic       
#    FASTFLAG = -O -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math                             
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
    PGLIBS   = -lgfortran -L/usr/local -lpgplot -L/usr/X11R6/lib -lX11 -L/scisoft/i386/lib -lpng12 -laquaterm /System/Library/Frameworks/Cocoa.framework/Cocoa
#   PGLIBS   = -L/usr/local/lib -lpgplot -L/usr/X11R6/lib -lX11 -lpng12 -laquaterm /System/Library/Frameworks/Cocoa.framework/Cocoa
#   PGLIBS   = -L/usr/local/lib /usr/local/lib/libpgplot.a -L/usr/X11R6/lib -lX11 -lpng12 -laquaterm /System/Library/Frameworks/Cocoa.framework/Cocoa
    CFITSIOLIBS = -L/star/lib/ -lcfitsio #-lmkl_core -lmkl_intel
#    PGLIBS   = -L/sw/lib/pgplot -lpgplot -L/usr/X11R6/lib -lX11 -L/sw/lib -lpng -laquaterm /System/Library/Frameworks/Cocoa.framework/Cocoa
#    CFITSIOLIBS = -L/Users/${USER}/cfitsio/lib -lcfitsio

#    STARLIBS = -L/star/lib `ndf_link`                                                                   
#    F2CLIB   = -L$(PATH2G2C) -lg2c                                                                      
    EXE = torus.intelmac
    MK_MPI_SRC = no
    KNOWN_SYSTEM = yes
endif



ifeq ($(SYSTEM), g95)
    F95COMP  = g95
    F77COMP  = g95
    F77FLAG  = -w -fno-second-underscore
    DEBUGFLAG= -g -ftrace=full
    FASTFLAG = -O3 -mcpu=G5
#    FASTFLAG = -O -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math	
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    PGLIBS   = -L/star/lib `pgplot_link`
    CFITSIOLIBS = -lcfitsio
    STARLIBS = -L/star/lib `ndf_link`
#    F2CLIB   = -L$(PATH2G2C) -lg2c
    EXE = torus.g95
    MK_MPI_SRC = no
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
    PGLIBS   = -L/star/lib `pgplot_link`
    STARLIBS = -L/star/lib `ndf_link`
    F2CLIB   = -L$(PATH2G2C) -lg2c
    EXE = torus.sun
    MK_MPI_SRC = no
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
#    PGLIBS   = -L/star/lib `pgplot_link`
    PGLIBS   = -L/sw/lib/pgplot -lpgplot
#    STARLIBS = -L/star/lib `ndf_link` 
    STARLIBS = -L/usr/X11R6/lib -lX11
    F2CLIB   = -L/usr/local/lib -lg2c
    EXE = torus.osx
    MK_MPI_SRC = no
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
    PGLIBS   = -L/star/lib `pgplot_link`
    STARLIBS = -L/star/lib `ndf_link`
#    STARLIBS = -L/Users/th/osxpgplot -lpgplot -L/usr/X11R6/lib -lX11
    F2CLIB   = -L/usr/local/lib -lg2c
    EXE = torus.gfortran
    MK_MPI_SRC = no
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
    STARLIBS = /local-star/lib/libpgplot.a  -L/Users/th/macmpi -lmacmpi /System/Library/Frameworks/Carbon.framework/Carbon -L/usr/lib -lSystemStubs /usr/X11R6/lib/libX11.a /usr/local/lib/libg2c.a /star/lib/libndf.a /star/lib/libary.a /star/lib/libhds.a /star/lib/liberr_standalone.a /star/lib/libpsx.a /star/lib/libchr.a /star/lib/libprm.a /star/lib/libprm_a.a /star/lib/libast.a /star/lib/libast_wcslib.a /star/lib/libast_slalib.a /star/lib/libsla.a /star/lib/libast_grf_2.0.a /star/lib/libast_grf_3.2.a /star/lib/libast_pass2.a /star/lib/libast_ems.a /star/lib/libemsf.a /star/lib/libems.a /star/lib/libcnf.a 
#    STARLIBS = -L/Users/th/osxpgplot -lpgplot -L/usr/X11R6/lib -lX11
#    F2CLIB   = -L/usr/local/lib -lg2c
    EXE = torus.macmpi
    MK_MPI_SRC = yes
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
    PGLIBS   = -L/star/lib `pgplot_link`
    STARLIBS = -L/star/lib `ndf_link` 
    F2CLIB   = -L$(PATH2G2C) -lg2c
    EXE = torus.g95
    MK_MPI_SRC = no
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
    PGLIBS   = -L/star/lib `pgplot_link`
    F2CLIB   = -lg2c
    STATICLIBS = -lpthread
    EXE = torus.linux
    MK_MPI_SRC = no
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
    PGLIBS   = -L/star/lib `pgplot_link`
    STARLIBS = -L/star/lib `ndf_link`
    EXE = torus
    MK_MPI_SRC = no
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
    PGLIBS   = -L${GUANINE_HOME}/csr201/local-ppc/opt/pgplot -lpgplot
    CFITSIOLIBS = -L${GUANINE_HOME}/tjharrie/cfitsio -lcfitsio
    STARLIBS = -L/usr/X11R6/lib -lX11 -L/usr/lib -lSystemStubs 
    EXE = torus.lamosx
    MK_MPI_SRC = yes
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
    PGLIBS   = -L${GUANINE_HOME}/csr201/local-i386/opt/pgplot -lpgplot
    CFITSIOLIBS = -L${GUANINE_HOME}/tjharrie/cfitsio.i386 -lcfitsio
    STARLIBS = -L/usr/X11R6/lib -lX11 -L/usr/lib -lSystemStubs 
    EXE = torus.lamosx-intel
    MK_MPI_SRC = yes
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
    PGLIBS   = -L${GUANINE_HOME}/csr201/local-ppc/opt/pgplot -lpgplot
    CFITSIOLIBS = -L${GUANINE_HOME}/tjharrie/cfitsio -lcfitsio
    STARLIBS = -L/usr/X11R6/lib -lX11 -L/usr/lib -lSystemStubs 
    EXE = torus.ompiosx
    MK_MPI_SRC = yes
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), ompi)
    F77COMP  = mpif90
    F95COMP  = mpif90
    F95FLAG  = 
    F77FLAG  = 
    DEBUGFLAG= -g -ftrace=full -fbounds-check -fimplicit-none -Wall -Wno=140,141,112,155
    FASTFLAG = -O0
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.miracle
    SPECF77  = spec_f77.f.miracle
    PGLIBS   =   -L/sw/lib -lpgplot -lpng -laquaterm  /System/Library/Frameworks/Cocoa.framework/Cocoa
    CFITSIOLIBS =  -L/Users/acreman/cfitsio/lib -lcfitsio
    STARLIBS = -L/usr/X11R6/lib -lX11
    EXE = torus.ompi
    MK_MPI_SRC = yes
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
    PGLIBS   = -L${GUANINE_HOME}/csr201/local-i386/opt/pgplot -lpgplot
    CFITSIOLIBS = -L${GUANINE_HOME}/tjharrie/cfitsio.i386 -lcfitsio
    STARLIBS = -L/usr/X11R6/lib -lX11 -L/usr/lib -lSystemStubs 
    EXE = torus.ompiosx-intel
    MK_MPI_SRC = yes
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
    PGLIBS   = -L/star/lib `pgplot_link`
    STARLIBS = -L/star/lib `ndf_link` 
    F2CLIB   = ${MPI_INC} -L${TOP}/lib -lmpi -llamf77mpi -L$(PATH2G2C) -lg2c -lm
    EXE = torus.mpilam
    MK_MPI_SRC = yes
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
    PGLIBS   = -L/star/lib `pgplot_link`
    STARLIBS = -L/star/lib `ndf_link`
    F2CLIB   = ${MPI_INC} -L${TOP}/lib -lmpichf90 -lmpich -L$(PATH2G2C) -lg2c -lm
    EXE = torus.mpich
    MK_MPI_SRC = yes
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), linuxintel)
    F95COMP  = ifort 
    F77COMP  = ifort
    F77FLAG  = -w
    DEBUGFLAG= -g



#    FASTFLAG = -g -openmp -fpp -auto -I/usr/include/nptl -L/usr/lib/nptl
#    FASTFLAG = -I/usr/include/nptl -L/usr/lib/nptl -zero 
    UNIXMOD  = unix_mod.f90.linuxintel
    KINDMOD  = kind_mod.f90.linuxintel
    IMAGEF77 = image_f77.f.default
    SPECF77  = spec_f77.f.default
    PGLIBS   = -L/star/lib `pgplot_link`
    CFITSIOLIBS = -lcfitsio
    STARLIBS = -L/star/lib `ndf_link` 
    F2CLIB   = -L/usr/lib -lg2c
    EXE = torus.linuxintel
    MK_MPI_SRC = no
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), zen)

    F95COMP= ifort -I/sw/sdev/intel-mpi/v3.0.021/x86_64/include64
    F77COMP= ifort -I/sw/sdev/intel-mpi/v3.0.021/x86_64/include64
    FASTFLAG= -O3 -xT -mcmodel=medium -i-dynamic
    DEBUGFLAG=-g -O0 -mcmodel=medium -warn -warn nounused
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
    STARLIBS= -DMPI -L/sw/sdev/intel-mpi/v3.0.021/x86_64/lib64 -Xlinker -rpath -Xlinker $libdir -Xlinker -rpath -Xlinker /opt/intel/mpi-rt/3.0 -lmpi -lmpiif -lmpigi -lrt -lpthread -ldl
#    -L/star/lib `ndf_link` 
    PGLIBS= -L/home/tjharrie/libpgplot -lpgplot -lpng
#    -L/star/lib `pgplot_link`
    CFITSIOLIBS = -L/home/tjharrie/cfitsio -lcfitsio
    MK_MPI_SRC = yes
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
    PGLIBS   = -L/star/lib `pgplot_link`
    STARLIBS = -L/star/lib `ndf_link` 
    F2CLIB   = ${MPI_INC} -L${TOP}/lib -lmpi -llamf77mpi -L$(PATH2G2C) -lg2c -lm
    EXE = torus.mpiintel
    MK_MPI_SRC = yes
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
    PGLIBS   = -L/star/lib `pgplot_link`
    CFITSIOLIBS = -lcfitsio
    STARLIBS = -L/star/lib `ndf_link` 
    F2CLIB   = ${MPI_INC} -L$(PATH2G2C) -lg2c
    EXE = torus.openmpi
    MK_MPI_SRC = yes
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
    PGLIBS   = -L/home/rk/pgplot -lpgplot -lX11
    EXE = torus.miracle
    MK_MPI_SRC = no
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
    PGLIBS   = -L/home/th/pgplot -lpgplot
    CFITSIOLIBS = -L/home/th/cfitsio -lcfitsio
    EXE = torus.ukaff1a
    MK_MPI_SRC = yes
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
    PGLIBS   = -L/home/th/pgplot -lpgplot
    EXE = torus.omp
    MK_MPI_SRC = no
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
    PGLIBS   = -L/usr/lib -lxlf90 -lm -lc -lX11 -L${HOME}/software/pgplot -lpgplot
    EXE = torus.aix
    MK_MPI_SRC = no
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
    PGLIBS   = -L/star/lib `pgplot_link`
    CFITSIOLIBS = -lcfitsio
    STARLIBS = -L/star/lib `ndf_link` 
    F2CLIB   = ${MPI_INC} -L${TOP}/lib -lmpi -llamf77mpi -L$(PATH2G2C) -lg2c -lm
    EXE = torus.mpig95
    MK_MPI_SRC = yes
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
    PGLIBS   = -L$(HOME)/pgplot -lpgplot
    STARLIBS = -L/usr/X11R6/lib -lX11 -L/usr/lib -lSystemStubs
    EXE = torus.darwin
    MK_MPI_SRC = no
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
    PGLIBS   = -L${HOME}/pgplot -lpgplot
    STARLIBS = -L/usr/X11R6/lib -lX11 -L/usr/lib -lSystemStubs 
    EXE = latticeQED
    MK_MPI_SRC = yes
    MPIFLAG = -DMPI
    KNOWN_SYSTEM = yes
endif

#
# Deal with make options
#
ETCFLAG  =

ifeq ($(pgplot), no)
   PGSTUBS = pgplot_stubs.o
   PGLIBS  = 
else 
   PGSTUBS = 
endif

ifeq ($(cfitsio), no)
   CFITSIOSTUBS = cfitsio_stubs.o
   CFITSIOLIBS  =
else 
   CFITSIOSTUBS =
endif

ifeq ($(debug), yes)
   FASTFLAG =
   ETCFLAG += $(DEBUGFLAG)
endif

ifeq ($(profile), yes)
   ETCFLAG += -pg
endif

ifeq ($(static), yes)
   ETCFLAG += $(STATICFLAG)
else
   STATICLIBS =
endif

ifeq ($(sph), yes)
   SPHFLAG = -DSPH
else
   SPHFLAG = 
endif

#
# Source files
#
F90SOURCE = kind_mod.f90 constants_mod.f90 atom_mod.f90 vector_mod.f90  \
            unix_mod.f90 inputs_mod.f90 utils_mod.f90 \
            gaussian_mod.f90 grid_mod.f90 math_mod.f90 \
            phasematrix_mod.f90 photon_mod.f90 blob_mod.f90 \
            distortion_mod.f90 image_mod.f90 TTauri_mod.f90 puls_mod.f90 \
            disc_hydro_mod.f90 gas_opacity_mod.f90 \
            integratePath.f90 readIntPro.f90 mieDistPhaseMatrix.f90 \
            gridtype_mod.f90 getRefractiveIndex.f90 draineCrossSection.f90 \
            readTioCrossSection.f90 rayleighCrossSection.f90 \
            fillGridTio.f90 fillGridRayleigh.f90 fillGridMie.f90 \
            intersectCube.f90 fillGridThomson.f90 cont_read.f90 \
            getMeanMass.f90 input_variables.f90 amr_mod.f90 octal_mod.f90 \
            parameters_mod.f90 jets_mod.f90 \
            dust_mod.f90 source_mod.f90 spectrum_mod.f90 ion_mod.f90 \
            wr104_mod.f90 opacity_lte_mod.f90 density_mod.f90 \
            pixplot_module.f90 sph_data_class.f90 cluster_class.f90 \
            linked_list_class.f90 timing.f90  isochrone_class.f90 \
            filter_set_class.f90 ostar_mod.f90  \
            cluster_utils.f90 surface_mod.f90  math_mod2.f90 disc_class.f90 \
            hyd_col_coeff.f90 rotation_variables.f90 function_variables.f90 \
            flow_speed_variables.f90 clump_mod.f90 \
            discwind_class.f90 zeus_data_class.f90 luc_cir3d_class.f90 \
            jet_class.f90 cmfgen_class.f90  irreg_grid_data_class.f90 \
            modelatom_mod.f90  datacube_mod.f90 magField.f90 \
            romanova_class.f90 messages_mod.f90 starburst_mod.f90 \
            mpi_global_mod.f90 nrtype.f90 nrutil.f90

F77SOURCE = bhmie.f phfit2.f mieDistCrossSection.f palette.f \
            spec_f77.f image_f77.f 

PGSOURCE  = pgplot_stubs.f90

CAPF90SOURCE =	 photoionAMR_mod.F90 photoion_mod.F90 mpi_amr_mod.F90 hydrodynamics_mod.F90 cmf_mod.F90 \
	         diffusion_mod.F90 formal_solutions.F90 lucy_mod.F90 molecular_mod.F90 \
	         parallel_mod.F90 stateq_mod.F90 torusMain.F90

# Object files
#
F90OBJ = $(F90SOURCE:.f90=.o)
CAPF90OBJ =  $(CAPF90SOURCE:.F90=.o)
F77OBJ = $(F77SOURCE:.f=.o)
PGOBJ  = $(PGSOURCE:.f90=.o)
ALLOBJ = $(CAPF90OBJ) $(F90OBJ)  $(F77OBJ) $(PGOBJ)

#
# Compilers and compiler flags
#
F90 = $(F95COMP) $(F95FLAG) $(ETCFLAG) $(FASTFLAG) $(MATHTRAP) $(MPI_INC) $(MPIFLAG) $(SPHFLAG)
F77 = $(F77COMP) $(F77FLAG) $(ETCFLAG) $(FASTFLAG) $(MATHTRAP) $(MPI_INC)

#
# Special targets
#   - New .f90.raw files need to be added to both preprocess *and* distclean.
#
all: depends $(ALLOBJ)
    ifeq ($(KNOWN_SYSTEM), yes)
	    @echo ""
	    @echo "Linking object files and libraries..."
	    $(F90) -o $(EXE) \
	    $(F77OBJ) $(F90OBJ) $(CAPF90OBJ) \
	    $(PGSTUBS) $(PGLIBS) $(CFITSIOSTUBS) $(CFITSIOLIBS) $(STARLIBS) $(F2CLIB) $(STATICLIBS)
	    @echo "...done"
    else
	    @echo ""
	    @echo "make: WARNING: value of SYSTEM = $(SYSTEM) not recognised..."
	    @echo "               ...setting SOURCE_VARIANT to unknown"
    endif

lib: depends $(ALLOBJ)
    ifeq ($(KNOWN_SYSTEM), yes)
	    @echo ""
	    @echo "Linking object files making a library..."
	    $(AR) $(ARFLAGS) libtorus.a $(F90OBJ) $(CAPF90OBJ) $(F77OBJ)
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
	@ctags --langmap='fortran:+(*.f90.raw)' \
	       *.f90 *.f *.f90.raw

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
