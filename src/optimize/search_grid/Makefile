FC = gfortran
#FFLAGS = -fopenmp -O2 -ffast-math -w
FFLAGS = -O2 -ffast-math -w
FILES_FS = mopti_fs.f90
#FILES_SC = mopti_sc.f90
#FILES_NF = mopti_nf.f90

all: clean

mopti_fs: $(FILES_FS)
	$(FC) $(FFLAGS) -o $@ $(FILES_FS)

#mopti_sc: $(FILES_SC)
#	$(FC) $(FFLAGS) -o $@ $(FILES_SC)

#mopti_nf: $(FILES_NF)
#	$(FC) $(FFLAGS) -o $@ $(FILES_NF)

#clean: mopti_fs mopti_sc mopti_nf
clean: mopti_fs
	#rm *.mod
