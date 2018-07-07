############   MAKEFILE   #############

# set compiler
FC = gfortran

# command file
CMD = growclust

# compiler flags for f90 codes
#  (use -O for optimization)
#  (use -O -mcmodel=medium or -O -mcmodel=large for large datasets)
#  (use -fcheck=bounds for testing)
CFLAGS = -O
#CFLAGS = -O -mcmodel=medium
#CFLAGS = -Wall -fcheck=bounds

# object files to link to command file
OBJS  =  input_subs.o \
         vel_subs.o \
         stats_subs.o
         
# modules
MODS  = grow_params.mod

#########################################

# make all
all: clean $(MODS) $(OBJS) $(CMD)

# make clean
clean:
	\rm -f $(CMD) *.o *.mod
	
#########################################

# make command file
#	 (note $@ selects target, $(CMD) )
$(CMD): $(CMD).f90
	$(FC) $(CFLAGS) $(OBJS) $< -o $@
	
# make all modules
# 	(note $< selects first prerequisite, %.f90. note also that $(@F:.o=.f90) would work...)	
%.mod: %.f90
	$(FC) -c $(CFLAGS) $<
	
# make all .o files (%.o: %.f90)	
%.o: %.f90
	$(FC) $(CFLAGS) -c $< -o $@
