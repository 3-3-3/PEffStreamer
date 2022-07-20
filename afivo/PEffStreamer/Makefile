AF_DIR  := ..
LIBDIRS := $(AF_DIR)/external_libraries/silo/lib $(AF_DIR)/external_libraries/hypre/lib
LIBS	:= afivo silo HYPRE
OUTDIR	:= output

include $(AF_DIR)/src/makerules.make

.PHONY: all clean run run_2d run_3d

PROG := streamer 

all:	$(PROG)

clean:
	$(RM) $(PROG) *.o *.mod

run:	$(PROG)
	do\
		echo "Running $$prog";\
		./$$prog || { echo "$$prog failed" ; exit 1; };\
	done

$(OUTDIR):
	mkdir -p $@

# Set flags
$(PROG): $(AF_DIR)/lib_2d/libafivo.a | $(OUTDIR)
$(PROG): INCDIRS+=$(AF_DIR)/lib_2d
$(PROG): LIBDIRS+=$(AF_DIR)/lib_2d
$(PROG): FFLAGS += -DNDIM=2

# How to create executables
%_2d: %.f90
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
