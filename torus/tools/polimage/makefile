
COBJ =  polimage.o

EXE=	./polimage
#	F90=	f90 -w=all -O3 -ieee=full

F90=	f90 -g
F77=	f77 -g -w -fno-second-underscore


.f.o:
	$(F77) -c $< -o $(<:.f=.o)

LIBDIR= -L/star/lib
LIBS=	`pgplot_link` `ndf_link` -last -last_err -last_ems \
	-last_grf -last_pass2 \
	-last_pgplot -last_wcslib -last_slalib -lsla

all: polimage polimage2

polimage:	$(COBJ)
	$(F77)  -o $@ \
	$(COBJ) $(LIBDIR) $(LIBS)  -lg2c

polimage2:	polimage2.o
	$(F77)  -o $@ \
	polimage2.o  $(LIBDIR) $(LIBS)  -lg2c

clean:
	/bin/rm *.o *~ *.mod polimage polimage2


