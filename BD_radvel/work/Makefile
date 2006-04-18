F90 = g95 -g 
F77 = f77 -g -w

OBJ = chi_square.o 

EXE = rvxs

LIBDIR= 
LIBS=   
all: rvxs

rvxs: rv_chisq.o $(OBJ)
	$(F90)   rv_chisq.o  $(OBJ)  $(LIBDIR) $(LIBS) -o $(EXE)

rv_chisq.o : chi_square.o 

.%.o : %.f
	$(F77) -c -fno-second-underscore $<

%.o : %.f90
	$(F90) -c $< 

read_ndf.o: read_ndf.f
	$(F77) -c -fno-second-underscore $<


clean:
	rm -f *.o *.mod
