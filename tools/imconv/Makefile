COBJ =  imconv.o subs.o

EXE=	./imconv
#	F90=	f90 -w=all -O3 -ieee=full

F90=	g95 -g 
F77=	g95 -g -w -fno-second-underscore



LIBDIR= -L/star/lib
LIBS=	`pgplot_link` `ndf_link` -last -last_err -last_ems \
	-last_pass2 \
	-last_pgplot -last_wcslib -last_slalib -lsla

imconv: $(COBJ)
	$(F90)  -o $@ \
	$(COBJ) $(LIBDIR) $(LIBS)  -lg2c


# General rule for compilations
%.o : %.f90
	$(F90) -c  $<
%.o : %.f
	$(F77) -c  $<




clean:
	/bin/rm *.o *~ *.mod


