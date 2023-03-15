#LDFLAGS += -L$(hdf5_dir)/lib -L$(paramesh_dir)/libs -L$(phase_dir)/libs 
LDLIBS += -lgfortran -lparamesh -lmpi_paramesh -lmodules -lhdf5 -lhdf5_hl 

#FFLAGS+=-g

LDFLAGS += $(ADD_LIB)
LDFLAGS += $(LDLIBS)
LDFLAGS += -lparamesh

# List all application source files
sources := amr_multigrid_user_modules.f90 \
amr_multigrid_user_edits_03.f90 amr_1blk_bcset.f90 amr_initial_soln.f90 \
amr_multigrid_general_03.f90 amr_multigrid_stat_03.f90 \
amr_multigrid_control_03.f90 amr_multigrid_test_refinement.f90 \
non_linear_mg.f90 ceg_mesh_partitioning.f90 ceg_shampine_conv.f90 \
ceg_local_timestep.f90
csources := profilingLinker.c
#amr_test_refinement.f90 amr_advance_soln.f90 amr_timestep.f90 
%.o:%.f90 amr_multigrid_user_modules.f90 Makefile.gnu
	$(FC) -c $(FFLAGS) $<
%.o:%.c Makefile.gnu
	$(CC) -c $(CFLAGS) $<
objects := $(sources:.f90=.o)
exobjects := $(csources:.c=.o)


# Identify the main program
mainobject := $(.f90=.o)
exobject := $(.c=.o)

# Set the executable name
CMD := ../pf.intel


# compiles the program.
$(CMD): $(mainobject) $(objects) $(exobjects)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)


.PHONY: clean
clean: 
	$(RM) $(CMD) *.o *.i *.mod *~
