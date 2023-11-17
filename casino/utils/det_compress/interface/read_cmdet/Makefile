# Simple Makefile for read demo
default:
	gfortran -O3 -c casl.f90 ;\
        gfortran -O3 -o read_demo read_demo.f90 casl.o

clean:
	rm -f casl.o casl.mod read_demo
