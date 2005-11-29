
OBJ =  subs.o plot_discwind.o

EXE=	pdw

F90=	ifort -w -g
F77=	ifort -w -g -fno-second-underscore

%.o : %.f90
	$(F90) -c  $<
%.o : %.f
	$(F77) -c  $<

pdw:  $(OBJ)
	$(F90)  -o $@ \
	$(OBJ) 

subs.o: subs.f90
plot_discwind.o: plot_discwind.f90
plot_discwind.o: subs.o

clean:
	/bin/rm *.o *~ *.mod 


