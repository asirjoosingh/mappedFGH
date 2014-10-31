
include Makefile.gnu

EXE = mappedfgh.exe

$(EXE) : cubicspline.o schroedinger.o primat.o mappedfgh.o
	$(FC90) $(OPTS) -o $(EXE) mappedfgh.o schroedinger.o primat.o cubicspline.o $(LIB)

%.o : %.f90
	$(FC90) $(OPTS) -c $<

%.mod : %.f90
	$(FC90) $(OPTS) -c $<

mappedfgh.o : cubicspline.mod mappedfgh.f90
	$(FC90) $(OPTS) -c mappedfgh.f90

primat.o : primat.f
	$(FC77) $(OPTS) -c primat.f

clean :
	rm *.o *.mod $(EXE)

