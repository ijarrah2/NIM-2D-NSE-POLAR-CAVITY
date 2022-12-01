# Fortran compiler
F90=gfortran
SILO_path=3rd_party/silo-4.11
# ============================================
F90SRC=main.f90
SRCOBJ=mesh.f90 user_parameters.f90 parameters.f90 general.f90 restart.f90 curvilinear.f90 geometry.f90 nim.f90 bc.f90 postprocess.f90 nse.f90
MODOBJ=$(SRCOBJ:.f90=.o) 
PROGRAM=nse2d
SILO_option=-I$(SILO_path)/include/ -L$(SILO_path)/lib/
OPTIONS= $(SILO_option) -lsilo -O3 -ffree-line-length-300 -fcheck=bounds -g -fopenmp -fdefault-real-8
# ============================================
$(PROGRAM):	$(MODOBJ)
	@echo $(MODOBJ)
	$(F90) $(OPTIONS) $(F90SRC) -lsilo -o $(PROGRAM)
%.o:	%.f90
	$(F90) -c $(OPTIONS) $^ -o $@
clean:	
	rm -f $(MODOBJ) $(PROGRAM) *.mod plots.visit plots/step*
cleanall:
	rm -f $(MODOBJ) $(PROGRAM) *.mod plots.visit plots/step*
cleansilo:
	rm -r 3rd_party/silo-4.11
	rm 3rd_party/silo-4.11.tar.gz
