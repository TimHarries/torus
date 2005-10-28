#----------------------------------------------------------------
# Usage :: 
#    1. make             : Make torus with optimization flag(s)
#    2. make clean       : Clean up directory
#    3. make debug=yes   : Make torus with debug option,
#    4. make profile=yes : Make torus with profiler option
#    5. make static=yes  : Make torus with static linking
#-----------------------------------------------------------------

KNOWN_SYSTEM=no
debug = no
static = no
profile= no
#PATH2G2C=/usr/lib/gcc-lib/i386-redhat-linux/3.2.2
PATH2G2C=/usr/lib/gcc/i586-mandrake-linux-gnu/3.4.1/

ifeq ($(SYSTEM), g95)
    F95COMP=g95
    F77COMP=g95
    FASTFLAG= -O3
#    FASTFLAG= -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    DEBUGFLAG=-g
    STATICFLAG=
    MATHTRAP=
    UNIXMOD=unix_mod.f90.linuxintel
    KINDMOD=kind_mod.f90.linuxintel
    EXE=torus.g95
    F77FLAG= -w -fno-second-underscore
    IMAGEF77=image_f77.f.default
    SPECF77=spec_f77.f.default
    F2CLIB=-L$(PATH2G2C) -lg2c
    STATICLIBS=
    STARLIBS=-L/star/lib `ndf_link` `pgplot_link`
    MK_MPI_SRC = no
    KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), osx)
    F95COMP=g95
    F77COMP=g95
    FASTFLAG= 
#    FASTFLAG= -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    DEBUGFLAG=-g
    STATICFLAG=
    MATHTRAP=
    UNIXMOD=unix_mod.f90.linuxintel
    KINDMOD=kind_mod.f90.linuxintel
    EXE=torus.osx
    F77FLAG= -w -fno-second-underscore
    IMAGEF77=image_f77.f.default
    SPECF77=spec_f77.f.default
#    F2CLIB=-L$(PATH2G2C) -lg2c
    F2CLIB= -L/usr/local/lib -lg2c
    STATICLIBS=
    STARLIBS=-L/star/lib `ndf_link` `pgplot_link`
#    STARLIBS=-L/Users/th/osxpgplot -lpgplot -L/usr/X11R6/lib -lX11
    MK_MPI_SRC = no
    KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), gfortran)
    F95COMP=gfortran
    F77COMP=gfortran
    FASTFLAG= -O3
#    FASTFLAG= -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    DEBUGFLAG=-g
    STATICFLAG=
    MATHTRAP=
    UNIXMOD=unix_mod.f90.linuxintel
    KINDMOD=kind_mod.f90.linuxintel
    EXE=torus.gfortran
    F77FLAG= -w -fno-second-underscore
    IMAGEF77=image_f77.f.default
    SPECF77=spec_f77.f.default
#    F2CLIB=-L$(PATH2G2C) -lg2c
    F2CLIB= -L/usr/local/lib -lg2c
    STATICLIBS=
    STARLIBS=-L/star/lib `ndf_link` `pgplot_link`
#    STARLIBS=-L/Users/th/osxpgplot -lpgplot -L/usr/X11R6/lib -lX11
    MK_MPI_SRC = no
    KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), macmpi)
    F95COMP=g95 -I/Users/th/macmpi -fno-second-underscore
    F77COMP=g95 -fno-second-underscore
    FASTFLAG=
#    FASTFLAG= -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -ffast-math
    DEBUGFLAG=-g
    STATICFLAG= -static
    MATHTRAP=
    UNIXMOD=unix_mod.f90.linuxintel
    KINDMOD=kind_mod.f90.linuxintel
    EXE=torus.macmpi
    F77FLAG= -w -fno-second-underscore
    IMAGEF77=image_f77.f.default
    SPECF77=spec_f77.f.default
#    F2CLIB=-L$(PATH2G2C) -lg2c
#    F2CLIB= -L/usr/local/lib -lg2c
    STATICLIBS=
    STARLIBS= /local-star/lib/libpgplot.a  -L/Users/th/macmpi -lmacmpi /System/Library/Frameworks/Carbon.framework/Carbon \
-L/usr/lib -lSystemStubs /usr/X11R6/lib/libX11.a /usr/local/lib/libg2c.a \
/star/lib/libndf.a /star/lib/libary.a \
/star/lib/libhds.a /star/lib/liberr_standalone.a /star/lib/libpsx.a /star/lib/libchr.a /star/lib/libprm.a \
/star/lib/libprm_a.a /star/lib/libast.a /star/lib/libast_wcslib.a /star/lib/libast_slalib.a /star/lib/libsla.a \
/star/lib/libast_grf_2.0.a /star/lib/libast_grf_3.2.a /star/lib/libast_pass2.a /star/lib/libast_ems.a /star/lib/libemsf.a \
/star/lib/libems.a /star/lib/libcnf.a 

#    STARLIBS=-L/Users/th/osxpgplot -lpgplot -L/usr/X11R6/lib -lX11
    MK_MPI_SRC = yes
    KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), gfort)
    F95COMP=gfc
    F77COMP=gfc
#    FASTFLAG= -O3
    FASTFLAG= -O3 -mtune=athlon-mp -march=athlon-mp -msse -mmmx -static
    DEBUGFLAG=-g 
    STATICFLAG=
    MATHTRAP=
    UNIXMOD=unix_mod.f90.linuxintel
    KINDMOD=kind_mod.f90.linuxintel
    EXE=torus.g95
    F77FLAG= -w -fno-second-underscore
    IMAGEF77=image_f77.f.default
    SPECF77=spec_f77.f.default
    F2CLIB=-L$(PATH2G2C) -lg2c
    STATICLIBS=
    STARLIBS=-L/star/lib `ndf_link` `pgplot_link`
    MK_MPI_SRC = no
    KNOWN_SYSTEM=yes
endif

ifeq  ($(SYSTEM), linux)
#   F95COMP=f95 -maxcontin=2500
   F95COMP=f95
   F77COMP=f77
   FASTFLAG=-O3 -w 
   DEBUGFLAG=-g -C=all
   STATICFLAG=-Bstatic
   MATHTRAP=-ieee=stop
   UNIXMOD=unix_mod.f90.linux
   KINDMOD=kind_mod.f90.linux
   EXE=torus.linux
   F77FLAG=-fno-second-underscore  -w
   IMAGEF77=image_f77.f.default
   SPECF77=spec_f77.f.default
   F2CLIB=-lg2c
   STATICLIBS=-lpthread
   STARLIBS=-L/star/lib `ndf_link` `pgplot_link`
   MK_MPI_SRC = no
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), alpha)
   F95COMP=f90
   F77COMP=f77
   FASTFLAG=-fast
   DEBUGFLAG=-g
   STATICFLAG=
   MATHTRAP=-fpe3
   UNIXMOD=unix_mod.f90.osf
   EXE=torus
   F77FLAG= 
   IMAGEF77=image_f77.f.default
   SPECF77=spec_f77.f.default
   KINDMOD=kind_mod.f90.osf
   F2CLIB= 
   STATICLIBS=
   STARLIBS=-L/star/lib `ndf_link` `pgplot_link`
   MK_MPI_SRC = no
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), lamosx)
   TOP=/usr/local
   TOP2=/usr/local
   F95COMP=${TOP}/bin/mpif77
   F77COMP=${TOP2}/bin/mpif77
   FASTFLAG=-O3 -I${TOP}/include  -dcfuns  -w 
   DEBUGFLAG=-g -I${TOP}/include  -dcfuns 
   STATICFLAG=
#   MATHTRAP=-ieee=stop
   MATHTRAP=
   UNIXMOD=unix_mod.f90.linuxintel
   KINDMOD=kind_mod.f90.linuxintel
   EXE=torus.lamosx
   F77FLAG=-fno-second-underscore -w
   IMAGEF77=image_f77.f.default
   SPECF77=spec_f77.f.default
   F2CLIB=-I${TOP}/include -L${TOP}/lib  -lmpi -llamf77mpi -L$(PATH2G2C) -lg2c -lm -L/usr/lib -lSystemStubs
   STATICLIBS=
   STARLIBS=-L/star/lib `ndf_link` `pgplot_link`
   MK_MPI_SRC=yes
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), mpilam)
   TOP=/h/rk/mpi_nag
   TOP2=/h/rk/mpi_g77
   F95COMP=${TOP}/bin/mpif95 -mismatch
   F77COMP=${TOP2}/bin/mpif77
   FASTFLAG=-O3 -I${TOP}/include  -dcfuns  -w 
   DEBUGFLAG=-g -I${TOP}/include  -dcfuns 
   STATICFLAG=
#   MATHTRAP=-ieee=stop
   MATHTRAP=
   UNIXMOD=unix_mod.f90.linux
   KINDMOD=kind_mod.f90.linux
   EXE=torus.mpilam
   F77FLAG=-fno-second-underscore -w
   IMAGEF77=image_f77.f.default
   SPECF77=spec_f77.f.default
   F2CLIB=-I${TOP}/include -L${TOP}/lib  -lmpi -llamf77mpi -L$(PATH2G2C) -lg2c -lm
   STATICLIBS=
   STARLIBS=-L/star/lib `ndf_link` `pgplot_link`
   MK_MPI_SRC=yes
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), mpich)
   TOP=/h/th/mpich
   TOP2=/h/th/mpich_g77
   F95COMP=${TOP}/bin/mpif90 -mismatch
   F77COMP=${TOP2}/bin/mpif77
   FASTFLAG=-O3 -I${TOP}/include -dcfuns -w
   DEBUGFLAG=-g -I${TOP}/include -dcfuns
   STATICFLAG=
#   MATHTRAP=-ieee=stop
   MATHTRAP=
   UNIXMOD=unix_mod.f90.linux
   KINDMOD=kind_mod.f90.linux
   EXE=torus.mpich
   F77FLAG=-fno-second-underscore -w
   IMAGEF77=image_f77.f.default
   SPECF77=spec_f77.f.default
   F2CLIB=-I${TOP}/include -L${TOP}/lib  -lmpichf90 -lmpich -L$(PATH2G2C) -lg2c -lm
   STATICLIBS=
   STARLIBS=-L/star/lib `ndf_link` `pgplot_link`
   MK_MPI_SRC=yes
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), linuxintel)
    F95COMP=ifort
    F77COMP=ifort
#    FASTFLAG= -openmp -fpp -auto -I/usr/include/nptl -L/usr/lib/nptl
    FASTFLAG= -I/usr/include/nptl -L/usr/lib/nptl -zero
    DEBUGFLAG=-g -O0 -zero
    STATICFLAG=
    MATHTRAP=-mp
    UNIXMOD=unix_mod.f90.linuxintel
    KINDMOD=kind_mod.f90.linuxintel
    EXE=torus.linuxintel
    F77FLAG=-w
    IMAGEF77=image_f77.f.default
    SPECF77=spec_f77.f.default
    F2CLIB=-L$(PATH2G2C) -lg2c
    STATICLIBS=
    STARLIBS=-L/star/lib `ndf_link` `pgplot_link`
    MK_MPI_SRC = no
    KNOWN_SYSTEM=yes
endif


ifeq ($(SYSTEM), mpiintel)
    TOP=$(HOME)/mpi_intel
    TOP2=$(HOME)/mpi_g77
    F95COMP=${TOP}/bin/mpif77 -I${TOP}/include -mp 
    F77COMP=${TOP2}/bin/mpif77 -I${TOP2}/include 
    FASTFLAG= -O2 -I/usr/include/nptl -L/usr/lib/nptl -zero
    #DEBUGFLAG=-g -O0 -CA -CU -CV -CS -CB  
    DEBUGFLAG=-g -O0 -CA -CV -CS -CB -zero
    STATICFLAG=
    MATHTRAP= 
    UNIXMOD=unix_mod.f90.linuxintel
    KINDMOD=kind_mod.f90.linuxintel
    EXE=torus.mpiintel
    F77FLAG=-fno-second-underscore -w
    IMAGEF77=image_f77.f.default
    SPECF77=spec_f77.f.default
    F2CLIB=-I${TOP}/include -L${TOP}/lib  -lmpi -llamf77mpi \
           -L$(PATH2G2C) -lg2c -lm
    STATICLIBS=
    STARLIBS=-L/star/lib `ndf_link` `pgplot_link` 
    MK_MPI_SRC =yes
    KNOWN_SYSTEM=yes
endif


ifeq ($(SYSTEM), mpig95)
   TOP=/h/rk/mpi_g95
   TOP2=
   F95COMP=${TOP}/bin/mpif77 -w
   F77COMP=${TOP}/bin/mpif77 -w
   FASTFLAG=-O3 -I${TOP}/include  
   DEBUGFLAG=-g -I${TOP}/include  
   STATICFLAG=
   MATHTRAP=-ieee=stop
   UNIXMOD=unix_mod.f90.linuxintel
   KINDMOD=kind_mod.f90.linuxintel
   EXE=torus.mpig95
   F77FLAG=-fno-second-underscore -w
   IMAGEF77=image_f77.f.default
   SPECF77=spec_f77.f.default
   F2CLIB=-I${TOP}/include -L${TOP}/lib  -lmpi -llamf77mpi -L$(PATH2G2C) -lg2c -lm
   STATICLIBS=
   STARLIBS=-L/star/lib `ndf_link` `pgplot_link`
   MK_MPI_SRC=yes
   KNOWN_SYSTEM=yes
endif



ifeq ($(SYSTEM), miracle)
   F95COMP=f90
   F77COMP=f77
   FASTFLAG=-O -OPT:Olimit=0 -mp -64
#      FASTFLAG=-OPT:Olimit=0 -O -mp
   DEBUGFLAG=-g
   STATICFLAG=
   MATHTRAP= 
   UNIXMOD=unix_mod.f90.miracle
   KINDMOD=kind_mod.f90.miracle
   EXE=torus.miracle
   F77FLAG= 
   IMAGEF77=image_f77.f.miracle
   SPECF77=spec_f77.f.miracle
   F2CLIB= 
   STATICLIBS=
   STARLIBS=-L/home/rk/pgplot -lpgplot -lX11
   MK_MPI_SRC = no
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), ukaff1a)
   F95COMP=mpif90 -qsuffix=f=f90
   F77COMP=mpif77
#   FASTFLAG=-O -OPT:Olimit=0 -mp -64
#      FASTFLAG=-OPT:Olimit=0 -O -mp
   DEBUGFLAG=-g
   STATICFLAG=
   MATHTRAP= 
   UNIXMOD=unix_mod.f90.miracle
   KINDMOD=kind_mod.f90.miracle
   EXE=torus.ukaff1a
   F77FLAG= 
   IMAGEF77=image_f77.f.miracle
   SPECF77=spec_f77.f.miracle
   F2CLIB= 
   STATICLIBS=
   STARLIBS=-L/home/th/pgplot -lpgplot
   MK_MPI_SRC = yes
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), darwin)
    F95COMP=xlf95 -qsuffix=f=f90
    F77COMP=f77
    FASTFLAG=-O2
    DEBUGFLAG=-g
    STATICFLAG=
    MATHTRAP= 
    UNIXMOD=unix_mod.f90.darwin
    KINDMOD=kind_mod.f90.darwin
    EXE=torus.darwin
    F77FLAG= 
    IMAGEF77=image_f77.f.darwin
    SPECF77=spec_f77.f.darwin
    F2CLIB=
    STATICLIBS=
    STARLIBS=-L/Users/mbate/Documents/ryuichi/pgplot -lpgplot 
    MK_MPI_SRC = no
    KNOWN_SYSTEM=yes
endif


ifeq ($(SYSTEM), ceres)
    F95COMP=xlf95 -bloadmap:outp -qmoddir=/lscratch -I/lscratch -qsuffix=f=f90 -qnolm -O3 -qarch=pwr3 -qmaxmem=-1
    F77COMP=xlf -qmoddir=/lscratch -I/lscratch -qsuffix=f=f -qnolm -O3 -qarch=pwr3 -qmaxmem=-1
    FASTFLAG=
    DEBUGFLAG=-g
    STATICFLAG=
    MATHTRAP= 
    UNIXMOD=unix_mod.f90.aix
    KINDMOD=kind_mod.f90.aix
    EXE=torus.aix
    F77FLAG= 
    IMAGEF77=image_f77.f.miracle
    SPECF77=spec_f77.f.miracle
    F2CLIB= 
    STATICLIBS=
    STARLIBS=-L/usr/lib -lxlf90 -lm -lc -lX11 -L${HOME}/software/pgplot -lpgplot
    MK_MPI_SRC = no
    KNOWN_SYSTEM=yes
endif


# Setting additional other options for special cases
ifeq ($(debug), yes)
   ETCFLAG =$(DEBUGFLAG)
   FASTFLAG =
else 
   DEBUGFLAG=
   ETCFLAG =$(DEBUGFLAG)
endif

ifeq ($(profile), yes)
   ETCFLAG = -pg
endif

ifneq ($(static), yes)
   STATICFLAG=
   STATICLIBS=
endif
#
#
#

.SUFFIXES:
.SUFFIXES: .f .o .f90


F90SOURCE= kind_mod.f90 constants_mod.f90 atom_mod.f90 vector_mod.f90  \
 unix_mod.f90  inputs_mod.f90 stateq_mod.f90 utils_mod.f90   gaussian_mod.f90  \
 grid_mod.f90 math_mod.f90  \
 phasematrix_mod.f90 photon_mod.f90 blob_mod.f90 distortion_mod.f90 \
 image_mod.f90 TTauri_mod.f90 puls_mod.f90 disc_hydro_mod.f90 gas_opacity_mod.f90 \
 torusMain.f90 integratePath.f90 readIntPro.f90 \
 mieDistPhaseMatrix.f90 gridtype_mod.f90 \
 getRefractiveIndex.f90 draineCrossSection.f90 \
 readTioCrossSection.f90 rayleighCrossSection.f90 \
 fillGridTio.f90 fillGridRayleigh.f90 fillGridMie.f90 \
 intersectCube.f90 fillGridThomson.f90 cont_read.f90 getMeanMass.f90 \
 input_variables.f90 amr_mod.f90 octal_mod.f90 parameters_mod.f90 diffusion_mod.f90 \
 lucy_mod.f90 jets_mod.f90 dust_mod.f90 source_mod.f90 spectrum_mod.f90 ion_mod.f90 photoion_mod.f90 \
 wr104_mod.f90 opacity_lte_mod.f90 density_mod.f90 pixplot_module.f90 \
 sph_data_class.f90 cluster_class.f90 linked_list_class.f90 timing.f90 \
 isochrone_class.f90 filter_set_class.f90 ostar_mod.f90 formal_solutions.f90 \
 cluster_utils.f90 surface_mod.f90  math_mod2.f90 disc_class.f90 \
 hyd_col_coeff.f90 rotation_variables.f90 function_variables.f90 \
 flow_speed_variables.f90 clump_mod.f90 parallel_mod.f90 \
 discwind_class.f90 zeus_data_class.f90 luc_cir3d_class.f90 \
 jet_class.f90 cmfgen_class.f90 messages_mod.f90


F77SOURCE= image_f77.f mieDistCrossSection.f palette.f\
           spec_f77.f bhmie.f phfit2.f

F90OBJ= $(F90SOURCE:.f90=.o)

F77OBJ= $(F77SOURCE:.f=.o)

ALLOBJ= $(F90OBJ) $(F77OBJ)

F90=    $(F95COMP)   $(ETCFLAG) $(STATICFLAG) $(FASTFLAG) $(MATHTRAP)
F77=    $(F77COMP)   $(ETCFLAG) $(STATICFLAG) $(FASTFLAG) $(MATHTRAP) $(F77FLAG)

LIBDIR=
LIBS=   



all: preprocess depends $(ALLOBJ)
   ifeq ($(KNOWN_SYSTEM), yes)
	@echo ""
	@echo " Linking object files and libraries..........."
	$(F90) -o $(EXE) \
	$(MODOBJ) $(F77OBJ) $(F90OBJ) $(STARLIBS) $(F2CLIB) $(STATICLIBS)
	@echo " ... done"
   else
	@echo ""
	@echo "make: WARNING: value of SYSTEM = $(SYSTEM) not recognised..."
	@echo "             ...setting SOURCE_VARIANT to unknown"
   endif



preprocess:
	@echo ""
	@echo "Creating TORUS for $(SYSTEM) system..........."
	@echo ""
	@echo "Creating links for system dependent files........"
	ln -sf $(UNIXMOD) unix_mod.f90
	ln -sf $(KINDMOD) kind_mod.f90
	ln -sf $(IMAGEF77) image_f77.f
	ln -sf $(SPECF77) spec_f77.f
	@echo " ... done"


depends: preprocess
	@echo ""
	@echo "Making dependency file .................."
	@touch -c make.depends
	@./makedepend
	@echo " ... done"
	@echo ""


mpif.h :
	@echo " " 
	@echo "Including mpif.h in the following file............."
	@echo " " 




# General rules to make object files
%.o: %.f90
	$(F90) -c $<

%.o: %.f
	$(F77) -c $<


# Specific rules to make object files
lucy_mod.o : lucy_mod.f90.raw amr_mod.o constants_mod.o grid_mod.o phasematrix_mod.o diffusion_mod.o \
             source_mod.o spectrum_mod.o timing.o vector_mod.o parallel_mod.o dust_mod.o messages_mod.o
        ifeq ($(MK_MPI_SRC), yes)
	   sed 's/\!$$MPI//g' lucy_mod.f90.raw   > lucy_mod.f90
        else
	   sed '/\!$$MPI/d'   lucy_mod.f90.raw   > lucy_mod.f90
        endif
	$(F90) -c lucy_mod.f90

disc_hydro_mod.o : disc_hydro_mod.f90.raw lucy_mod.o amr_mod.o constants_mod.o grid_mod.o phasematrix_mod.o \
             source_mod.o spectrum_mod.o timing.o vector_mod.o messages_mod.o
        ifeq ($(MK_MPI_SRC), yes)
	   sed 's/\!$$MPI//g' disc_hydro_mod.f90.raw   > disc_hydro_mod.f90
        else
	   sed '/\!$$MPI/d'   disc_hydro_mod.f90.raw   > disc_hydro_mod.f90
        endif
	$(F90) -c disc_hydro_mod.f90

stateq_mod.o: stateq_mod.f90.raw cmfgen_class.o constants_mod.o grid_mod.o gridtype_mod.o \
	hyd_col_coeff.o input_variables.o integratePath.o jets_mod.o \
	math_mod.o opacity_lte_mod.o parallel_mod.o utils_mod.o vector_mod.o  messages_mod.o
        ifeq ($(MK_MPI_SRC), yes)
	   sed 's/\!$$MPI//g' stateq_mod.f90.raw > stateq_mod.f90
        else
	   sed '/\!$$MPI/d'   stateq_mod.f90.raw > stateq_mod.f90 
        endif
	$(F90) -c stateq_mod.f90


parallel_mod.o: parallel_mod.f90.raw  
        ifeq ($(MK_MPI_SRC), yes)
	   sed 's/\!$$MPI//g' parallel_mod.f90.raw > parallel_mod.f90
        else
	   sed '/\!$$MPI/d'   parallel_mod.f90.raw > parallel_mod.f90
        endif
	$(F90) -c parallel_mod.f90


torusMain.o: torusMain.f90.raw TTauri_mod.o amr_mod.o blob_mod.o cluster_class.o \
	cluster_utils.o cmfgen_class.o constants_mod.o disc_class.o \
	disc_hydro_mod.o discwind_class.o distortion_mod.o dust_mod.o \
	filter_set_class.o gas_opacity_mod.o grid_mod.o gridtype_mod.o \
	image_mod.o input_variables.o inputs_mod.o integratePath.o \
	isochrone_class.o jet_class.o jets_mod.o kind_mod.o luc_cir3d_class.o \
	lucy_mod.o math_mod.o phasematrix_mod.o photon_mod.o puls_mod.o \
	source_mod.o spectrum_mod.o sph_data_class.o stateq_mod.o \
	surface_mod.o timing.o unix_mod.o utils_mod.o vector_mod.o \
	wr104_mod.o  messages_mod.o photoion_mod.o


        ifeq ($(MK_MPI_SRC), yes)
	   sed 's/\!$$MPI//g' torusMain.f90.raw  > torusMain.f90
        else
	   sed '/\!$$MPI/d' torusMain.f90.raw  > torusMain.f90
        endif
	$(F90) -c torusMain.f90


clean:
	@/bin/rm -f *.o *~ *.mod *.il *.pc *.d *.rnreg *.life *.cfg *.cse *.combine


# including the dependency file
include make.depends



# DO NOT DELETE
