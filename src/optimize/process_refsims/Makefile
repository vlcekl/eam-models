FC = gfortran
FFLAGS = -O2 -ffast-math -w

FILES = trajecta.f90

all: clean

trajecta: $(FILES)
	$(FC) $(FFLAGS) -o $@ $(FILES)

clean: trajecta
	#rm *.mod
