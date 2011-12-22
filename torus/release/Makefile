memcheck = no
openmp  = no 

ifeq ($(SYSTEM), nagfor)
    F95COMP  = nagfor
    F95FLAG  = -colour -kind=byte -fpp -dcfuns -f2003 -DUSEIEEEISNAN -DUSEUNIXENV
    DEBUGFLAG= -C=all -g -gline -info -nocheck_modtime
    FASTFLAG = -O4
    CFITSIOLIBS = -L/Users/${USER}/cfitsio_nagfor/lib -lcfitsio 
    EXE = torus.nagfor
    KNOWN_SYSTEM = yes
    memcheck = no
endif

ifeq ($(SYSTEM), g95)
    F95COMP  = g95
    DEBUGFLAG= -g3 -fbounds-check -ftrace=full -fimplicit-none -freal=nan -Wall -Wunused-internal-procs -Wunused-parameter -Wno=140,141,112 -Werror
    FASTFLAG = -O3
    CFITSIOLIBS = -lcfitsio
    EXE = torus.g95
    KNOWN_SYSTEM = yes
endif

ifeq ($(SYSTEM), gfortran)
    ifeq ($(mpi), yes)
        F95COMP = mpif90
        MPIFLAG = -DMPI
    else
        F95COMP = gfortran
    endif
    DEBUGFLAG= -ggdb -fbacktrace -fbounds-check -Wunderflow -Wall -Werror -pedantic-errors -std=f2003 -fall-intrinsics -ffpe-trap=invalid,zero,overflow #,underflow,precision,denormal
    FASTFLAG = -O3
    EXE = torus.gfortran
    KNOWN_SYSTEM = yes
    CFITSIOLIBS = -L/Users/${USER}/cfitsio/lib -lcfitsio
endif

ifeq ($(SYSTEM), ifort)
    F95COMP= ifort
    FASTFLAG= -O2
    DEBUGFLAG=-g -O0 -warn -traceback -fpe0 -fp-stack-check -gen-interfaces -check all -check noarg_temp_created
    STATICFLAG=
    EXE=torus.ifort
    F2CLIB=
    STATICLIBS= -Bstatic
    STARLIBS= 
    KNOWN_SYSTEM=yes
endif

#
# Deal with make options
#
ETCFLAG  =

ifeq ($(memcheck), yes)
   ETCFLAG += -DMEMCHECK
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

TORUSMAINSOURCE= torusMainV2.F90
INPUTSMODSOURCE= inputs_modV2.f90

#
# Source files
#

F90SOURCE = torus_version_mod.f90 kind_mod.f90 constants_mod.f90 atom_mod.f90 vector_mod.f90  \
            utils_mod.f90  bhmie_mod.f90 mieDistCrossSection_mod.f90 \
            gaussian_mod.f90 grid_mod.f90 math_mod.f90 dimensionality_mod.f90 \
            phasematrix_mod.f90 photon_mod.f90 blob_mod.f90 \
            distortion_mod.f90  TTauri_mod.f90 \
            disc_hydro_mod.f90 gas_opacity_mod.f90 \
            integratePath.f90 readIntPro.f90 mieDistPhaseMatrix.f90 \
            gridtype_mod.f90 \
            readTioCrossSection.f90  \
            fillGridTio.f90 fillGridRayleigh.f90 \
            fillGridThomson.f90  \
            amr_mod.f90 amr_utils_mod.f90  \
            jets_mod.f90 h21cm_mod.f90 \
            dust_mod.f90 source_mod.f90 spectrum_mod.f90 \
            density_mod.f90 vh1_mod.f90 \
            timing.f90  isochrone_class.f90 \
            filter_set_class.f90 ostar_mod.f90  \
            surface_mod.f90  math_mod2.f90 disc_class.f90 \
            hyd_col_coeff.f90  clump_mod.f90 \
            discwind_class.f90 zeus_data_class.f90 luc_cir3d_class.f90 \
            cmfgen_class.f90 magField.f90 \
            romanova_class.f90 messages_mod.f90 starburst_mod.f90 \
            magnetic_mod.f90 sed_mod.f90 fits_utils_mod.f90 \
            image_utils_mod.f90 interferometry_mod.f90 units_mod.f90
            


CAPF90SOURCE =	 octal_mod.F90 mpi_global_mod.F90 random_mod.F90 unix_mod.F90  \
	         mpi_amr_mod.F90  \
	         diffusion_mod.F90 formal_solutions.F90 lucy_mod.F90  \
	         parallel_mod.F90 stateq_mod.F90 image_mod.F90 datacube_mod.F90 vtk_mod.F90 \
                 phaseloop_mod.F90 gridio_mod.F90  \
	         setupamr_mod.F90 physics_mod.F90 outputs_mod.F90 memory_mod.F90 $(TORUSMAINSOURCE)

# Object files
#
INPUTSMODOBJ = $(INPUTSMODSOURCE:.f90=.o)
F90OBJ = $(F90SOURCE:.f90=.o)
CAPF90OBJ =  $(CAPF90SOURCE:.F90=.o)
CAPF90OBJLIB =  $(CAPF90SOURCELIB:.F90=.o)
ALLOBJ = $(CAPF90OBJ) $(F90OBJ)  $(INPUTSMODOBJ)

#
# Compilers and compiler flags
#

F90 = $(F95COMP) $(F95FLAG) $(FASTFLAG) $(ETCFLAG)

#
#

all: depends $(ALLOBJ)
    ifeq ($(KNOWN_SYSTEM), yes)
	    @echo ""
	    @echo "Linking object files and libraries..."
	    $(F90) $(TRACEFLAG) -o $(EXE) \
	    $(INPUTSMODOBJ) $(F90OBJ) $(CAPF90OBJ) \
	    $(CFITSIOLIBS) $(STARLIBS) $(F2CLIB) $(STATICLIBS) $(MKLIBS) $(XMLLIBS)
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
