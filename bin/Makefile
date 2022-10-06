#############################################################################
# If you have problems with this makefile, contact Romain.Teyssier@gmail.com
#############################################################################
# Compilation time parameters
NVECTOR = 32
NDIM = 3
NHILBERT = 1
CONDINIT = coeur
UNITS = coeur
NPRE = 8
NVAR = 5
NENER = 0
COMPILER = GNU
DEBUG = 0
MPI = 0
MPIF90 = mpif90
SOLVER = hydro
PATCH = 
EXEC = ramses
#############################################################################
GITBRANCH = $(shell git rev-parse --abbrev-ref HEAD)
GITHASH = $(shell git log --pretty=format:'%H' -n 1)
GITREPO = $(shell git config --get remote.origin.url)
BUILDDATE = $(shell date +"%D-%T")
DEFINES = -DNVECTOR=$(NVECTOR) -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DNENER=$(NENER) -DNVAR=$(NVAR) \
	  -DNHILBERT=$(NHILBERT) -DHYDRO -DGRAV -DCONDINIT${CONDINIT} -DUNITS${UNITS}
#############################################################################
# GNU compiler (gfortran)
ifeq ($(COMPILER),GNU)
   FFLAGS = -x f95-cpp-input $(DEFINES)
   ifeq ($(MPI),1)
      F90 = $(MPIF90)
   else
      F90 = gfortran
      FFLAGS += -DWITHOUTMPI
   endif
   F90 += -frecord-marker=4 -fbacktrace -ffree-line-length-none -g -fimplicit-none
   ifeq ($(DEBUG),1)
      F90 += -O0 -fbounds-check -Wuninitialized -Wall
      FFLAGS += -ffpe-trap=zero,underflow,overflow,invalid -finit-real=nan
   else
      F90 += -O3
   endif
endif
#############################################################################
# Intel compiler
ifeq ($(COMPILER),INTEL)
   FFLAGS = -cpp $(DEFINES)
   ifeq ($(MPI),1)
      F90 = $(MPIF90)
   else
      F90 = ifort
      FFLAGS += -DWITHOUTMPI
   endif
   F90 += -fp-model source
   ifeq ($(DEBUG),1)
      F90 += -warn all -O0 -g -traceback -check bounds
      FFLAGS += -fpe0 -ftrapuv -init=zero -init=snan -init=arrays
   else
      F90 += -O3
   endif
endif
#############################################################################
C = ifort
CFLAGS = -g -O3
#############################################################################
MOD = mod
#############################################################################
# MPI librairies
LIBMPI = 
#LIBMPI = -lfmpi -lmpi -lelan

LIBS = $(LIBMPI)
#############################################################################
# Sources directories are searched in this exact order
VPATH = $(PATCH):../$(SOLVER):../hydro:../pm:../poisson:../amr:../mdl1
#############################################################################

# Modules
MODOBJ = amr_parameters.o hydro_parameters.o oct_commons.o hash.o pm_parameters.o pm_commons.o hilbert.o domain_hilbert.o \
         hydro_commons.o amr_commons.o poisson_parameters.o poisson_commons.o \
	 mdl_commons.o mdl.o ramses_commons.o call_back.o cache_commons.o cache.o gadgetreadfile.o write_makefile.o write_patch.o \
	 write_gitinfo.o output_hydro.o output_poisson.o output_part.o output_amr.o movie.o nbors_utils.o nbors_utils_p.o hydro_flag.o \
	 multigrid_fine_coarse.o rho_fine.o marshal.o upload.o init_hydro.o init_part.o init_amr.o smooth.o flag_utils.o \
	 init_time.o read_params.o input_part_grafic.o input_part_ascii.o input_part_restart.o input_part.o \
	 update_time.o timer.o init_refine_basegrid.o load_balance.o refine_utils.o \
	 init_refine_adaptive.o init_refine_restart.o title.o units.o move_fine.o rho_ana.o \
	 grav_ana.o phi_fine_cg.o interpol_phi.o force_fine.o poisson_flag.o multigrid_fine_commons.o \
	 input_hydro_condinit.o input_hydro_grafic.o init_flow_fine.o interpol_hydro.o \
	 condinit.o godunov_utils.o write_screen.o courant_fine.o newdt_fine.o \
	 godunov_fine.o umuscl.o uplmde.o synchro_hydro_fine.o cooling_fine.o \
	 task_manager.o amr_step.o adaptive_loop.o 

#############################################################################
ramses:	$(MODOBJ) ramses.o
	$(F90) $(MODOBJ) ramses.o -o $(EXEC)$(NDIM)d $(LIBS)
	rm write_makefile.f90
	rm write_patch.f90

amr2map:   amr_parameters.o timer.o hilbert.o domain_hilbert.o amr_commons.o hydro_parameters.o hydro_commons.o ../utils/f90/amr2map.f90
	 $(F90) amr_parameters.o timer.o hilbert.o domain_hilbert.o amr_commons.o hydro_parameters.o hydro_commons.o ../utils/f90/amr2map.f90 -o amr2map

#############################################################################
write_gitinfo.o: FORCE
	$(F90) $(FFLAGS) -DPATCH=\'$(PATCH)\' -DGITBRANCH=\'$(GITBRANCH)\' -DGITHASH=\'"$(GITHASH)"\' \
-DGITREPO=\'$(GITREPO)\' -DBUILDDATE=\'"$(BUILDDATE)"\' -c ../amr/write_gitinfo.f90 -o $@		
write_makefile.o: FORCE
	../utils/scripts/cr_write_makefile.sh $(MAKEFILE_LIST)
	$(F90) $(FFLAGS) -c write_makefile.f90 -o $@
write_patch.o: FORCE
	../utils/scripts/cr_write_patch.sh $(PATCH)
	$(F90) $(FFLAGS) -c write_patch.f90 -o $@
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@
%.o:%.c
	$(C) $(CFLAGS) -c $^ -o $@

FORCE:
#############################################################################
clean :
	rm *.o *.$(MOD)
#############################################################################
