.SUFFIXES:
.SUFFIXES: .f .o .f90

F90SOURCE= kind_mod.f90 constants_mod.f90 atom_mod.f90 vector_mod.f90  \
 unix_mod.f90 inputs_mod.f90 stateq_mod.f90 \
 utils_mod.f90 grid_mod.f90 math_mod.f90  \
 phasematrix_mod.f90 photon_mod.f90 blob_mod.f90 distortion_mod.f90 \
 image_mod.f90 TTauri_mod.f90 puls_mod.f90 torusMain.f90 \
 integratePath.f90   readIntPro.f90 \
 mieDistPhaseMatrix.f90 gridtype_mod.f90 \
 getRefractiveIndex.f90 draineCrossSection.f90 \
 readTioCrossSection.f90 rayleighCrossSection.f90 \
 fillGridTio.f90 fillGridRayleigh.f90 fillGridMie.f90 \
 intersectCube.f90 fillGridThomson.f90 cont_read.f90 getMeanMass.f90 \
 input_variables.f90 amr_mod.f90 octal_mod.f90 parameters_mod.f90 \
 lucy_mod.f90 jets_mod.f90 dust_mod.f90 source_mod.f90 spectrum_mod.f90 \
 wr104_mod.f90 opacity_lte_mod.f90 density_mod.f90 pixplot_module.f90 \
 sph_data_class.f90 cluster_class.f90 linked_list_class.f90 timing.f90 \
 isochrone_class.f90


F77SOURCE= image_f77.f mieDistCrossSection.f palette.f\
           spec_f77.f bhmie.f

F90OBJ= $(F90SOURCE:.f90=.o)

F77OBJ= $(F77SOURCE:.f=.o)

ALLOBJ= $(F90OBJ) $(F77OBJ)



F90=    $(F95COMP)   $(DEBUGFLAG) $(FASTFLAG) $(MATHTRAP)
F77=    $(F77COMP)   $(DEBUGFLAG) $(FASTFLAG) $(MATHTRAP) $(F77FLAG)

LIBDIR=
LIBS=   

.f90.o:
	$(F90) -c $<

.f.o:
	$(F77) -c $<


all:    $(ALLOBJ)
	$(F90) -o $(EXE) \
	$(MODOBJ) $(F77OBJ) $(F90OBJ) $(STARLIBS) $(F2CLIB) 



depends: 
	@./makedepend


clean:
	@/bin/rm -f *.o *~ *.mod *.il *.pc *.d ./make.depends


include make.depends

# DO NOT DELETE








