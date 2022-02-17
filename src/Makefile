all : enac

F90FLAGS = -cpp -g -fcheck=all
#FC       = ifort
FC       = gfortran

#NCINC    = -I/opt/netcdf/netcsdf-4.3.3.1/include/ -I/opt/netcdf/netcdf-fortran-4.4.2/include/ -I/opt/hdf5/hdf5-1.8.14/include/
#NCLIB    = -L/opt/netcdf/netcdf-4.3.3.1/lib/ -L/opt/netcdf/netcdf-fortran-4.4.2/lib/ -L/opt/hdf5/hdf5-1.8.14/lib/ -lnetcdf -lnetcdff -lhdf5

#NCINC    = -I/mingw64/include
#NCLIB    = -L/mingw64/lib/ -lnetcdf -lnetcdff -lhdf5


.SUFFIXES: .F90 .o
.F90.o:  
	$(FC) $(F90FLAGS) $(NCINC) $(NCLIB) -c $*.F90

.SUFFIXES: .f90 .o
.f90.o:  
	$(FC) $(F90FLAGS) $(NCINC) $(NCLIB) -c $*.f90

MSRCS = enac_commons.o enac_functions.o enac_inout.o enac_get_calcs.o enac_links.o enac_superindividuals.o enac_update_pop.o enac_decision_model.o enac_main.o

HSRCS =

FSRCS = 

SSRCS =

enac: $(MSRCS) $(HSRCS) $(FSRCS) $(SSRCS)
	$(FC) $(F90FLAGS) $(MSRCS) $(FSRCS) $(SSRCS) $(NCLIB) $(NCINC) -o enac.unix

clean:
	rm -f *.o *.mod *.unix core

