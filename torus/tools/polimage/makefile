
COBJ =  polimage.o

EXE=	./polimage
#	F90=	f90 -w=all -O3 -ieee=full

F90=	f90 -g
F77=	f77 -g -fno-second-underscore


.f.o:
	$(F77) -c $< -o $(<:.f=.o)

LIBDIR= -L/star/lib
LIBS=	`pgplot_link` `ndf_link` -last -last_err -last_ems \
	-last_grf -last_pass2 \
	-last_pgplot -last_wcslib -last_slalib -lsla

polimage:	$(COBJ)
	$(F77)  -o $@ \
	$(COBJ) $(LIBDIR) $(LIBS)  -lg2c

clean:
	/bin/rm *.o *~ *.mod polimage


