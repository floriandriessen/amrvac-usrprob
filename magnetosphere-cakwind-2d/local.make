# This is a sub-makefile for the main AMRVAC makefile
# All variables not defined in this local.make are inherited from makefile in
# the user problem folder

# Directory with extra modules for use in mod_usr.t
USR_SHAREDIR := ../share

vpath %.f90 $(USR_SHAREDIR)

# Local user fortran files
FSRCS := mod_kind_parameter.f90 mod_init_cakwind.f90
FOBJS := $(FSRCS:.f90=.o)
FMODS := $(FOBJS:.o=.mod)

$(info Extra user source file(s) included: $(FSRCS))

include $(USR_SHAREDIR)/share.dependencies

# Rule for the target
amrvac mod_usr.o: $(FOBJS)

# Build rule for local user object files
# Note amrvac's INC_DIRS because some user modules may depend on amrvac modules
%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS))

mod_%.mod: mod_%.f90 mod_%.o
	@test -f $@ || $(F90) $(F90FLAGS) -c $(@:.mod=.f90) -o $(@:.mod=.o) \
	$(addprefix -I,$(INC_DIRS))

clean_user:
	echo 'Cleaning local user objects'
	$(RM) $(FOBJS) $(FMODS)
