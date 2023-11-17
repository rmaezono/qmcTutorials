FC=ifort -O2
LIBDIR=$(ATSP)/lib
LIB=-lang
EXEC=${HOME}/bin/extractdet

OBJ = main.o extractcsf.o subshell.o error.o anglib.o

main: $(OBJ)
	 $(FC) -o $(EXEC)  $(OBJ)  $(LIB) -L$(LIBDIR)
# Here are the compile steps

main.o:main.f90 extractcsf.o error.o
	 $(FC) -c  main.f90

extractcsf.o:extractcsf.f90 subshell.o error.o
	 $(FC) -c  extractcsf.f90

subshell.o:subshell.f90 error.o anglib.o
	 $(FC) -c  subshell.f90

error.o:error.f90
	 $(FC) -c  error.f90

anglib.o:anglib.f90
	 $(FC) -c  anglib.f90


clean:
	 rm -f *.o *.mod
