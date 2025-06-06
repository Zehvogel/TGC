# WHIZARD: Makefile for process library 'default_lib'
# Automatically generated file, do not edit

# Integrity check (don't modify the following line!)
MD5SUM = '906D7A57E3E57789066650303B7E31D6'

# Library name
BASE = default_lib

# Compiler
FC = /cvmfs/sw.hsf.org/contrib/x86_64-almalinux9-gcc11.4.1-opt/gcc/14.2.0-yuyjov/bin/gfortran
CC = /cvmfs/sw.hsf.org/contrib/x86_64-almalinux9-gcc11.4.1-opt/gcc/14.2.0-yuyjov/bin/gcc

# Included libraries
FCINCL = -I/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/whizard/3.1.4-zjmc7r/lib/mod/whizard -I/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/whizard/3.1.4-zjmc7r/lib/mod/omega -I/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/whizard/3.1.4-zjmc7r/lib/mod/models -I/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/openloops/2.1.2-eh75iu/lib_src/openloops/mod 

# Compiler flags
FCFLAGS =  -fopenmp  -g -O2
FCFLAGS_PIC =  -fPIC
CFLAGS = -g -O2
CFLAGS_PIC = 
LDFLAGS = -L/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/whizard/3.1.4-zjmc7r/lib -lwhizard_main -lwhizard -lomega -I/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/libtirpc/1.3.3-db2ekl/include/tirpc -L/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/libtirpc/1.3.3-db2ekl/lib -ltirpc -Wl,-rpath,/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/hepmc3/3.3.0-hhj326/lib -L/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/hepmc3/3.3.0-hhj326/lib -lHepMC3  -Wl,-rpath,/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/lcio/2.22.6-zicklo/lib -Wl,-rpath,/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/lcio/2.22.6-zicklo/lib64 -L/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/lcio/2.22.6-zicklo/lib -L/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/lcio/2.22.6-zicklo/lib64 -llcio   -Wl,-rpath,/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/openloops/2.1.2-eh75iu/lib -L/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/openloops/2.1.2-eh75iu/lib -Wl,-rpath,/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/openloops/2.1.2-eh75iu/lib64 -L/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/openloops/2.1.2-eh75iu/lib64 -lopenloops    -lomega -L/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/whizard/3.1.4-zjmc7r/lib/whizard/models -lwhizard

# LaTeX setup
LATEX = no -halt-on-error
MPOST = no  -halt-on-error
DVIPS = no
PS2PDF = no
TEX_FLAGS = "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/whizard/3.1.4-zjmc7r/share/texmf/whizard:$$TEXINPUTS"
MP_FLAGS  = "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/whizard/3.1.4-zjmc7r/share/texmf/whizard:$$MPINPUTS"

# Libtool
LIBTOOL = /cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/whizard/3.1.4-zjmc7r/lib/whizard/libtool
FCOMPILE = @$(LIBTOOL) --silent --tag=FC --mode=compile
CCOMPILE = @$(LIBTOOL) --silent --tag=CC --mode=compile
LINK = @$(LIBTOOL) --silent --tag=FC --mode=link

# Compile commands (default)
LTFCOMPILE = $(FCOMPILE) $(FC) -c $(FCINCL) $(FCFLAGS) $(FCFLAGS_PIC)
LTCCOMPILE = $(CCOMPILE) $(CC) -c $(CFLAGS) $(CFLAGS_PIC)

# Default target
all: link diags

# Matrix-element code files
SOURCES += ww_i1.f90
OBJECTS += ww_i1.lo
ww_i1.f90:
	@echo  "  OMEGA     ww_i1.f90"
	@/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/whizard/3.1.4-zjmc7r/bin/omega_SM_ac.opt -o ww_i1.f90 -target:whizard -target:parameter_module parameters_SM_ac -target:module opr_ww_i1 -target:md5sum '84CE608809668E24C61FC5BB8254BA36' -target:openmp -fusion:progress -scatter 'e- e+ -> e-:e+ nue:nuebar u:ubar d:dbar' 
clean-ww_i1:
	@echo  "  RM        ww_i1.f90,.mod,.lo"
	@rm -f ww_i1.f90
	@rm -f opr_ww_i1.mod
	@rm -f ww_i1.lo
CLEAN_SOURCES += ww_i1.f90
CLEAN_OBJECTS += opr_ww_i1.mod
CLEAN_OBJECTS += ww_i1.lo
ww_i1.lo: ww_i1.f90
	$(LTFCOMPILE) $<

# Library driver
$(BASE).lo: $(BASE).f90 $(OBJECTS)
	$(LTFCOMPILE) $<
	@echo  "  FC       " $@

# Library
$(BASE).la: $(BASE).lo $(OBJECTS)
	@echo  "  FCLD     " $@
	$(LINK) $(FC) -module -rpath /dev/null $(FCFLAGS) $(LDFLAGS) -o $(BASE).la $^

# Main targets
link: compile $(BASE).la
compile: source $(OBJECTS) $(TEX_OBJECTS) $(BASE).lo
compile_tex: $(TEX_OBJECTS)
source: $(SOURCES) $(BASE).f90 $(TEX_SOURCES)
.PHONY: link diags compile compile_tex source

# Specific cleanup targets
clean-ww_i1:
.PHONY: clean-ww_i1

# Generic cleanup targets
clean-library:
	@echo  "  RM        $(BASE).la"
	@rm -f $(BASE).la
clean-objects:
	@echo  "  RM        $(BASE).lo $(BASE)_driver.mod $(CLEAN_OBJECTS)"
	@rm -f $(BASE).lo $(BASE)_driver.mod $(CLEAN_OBJECTS)
clean-source:
	@echo  "  RM        $(CLEAN_SOURCES)"
	@rm -f $(CLEAN_SOURCES)
clean-driver:
	@echo  "  RM        $(BASE).f90"
	@rm -f $(BASE).f90
clean-makefile:
	@echo  "  RM        $(BASE).makefile"
	@rm -f $(BASE).makefile
.PHONY: clean-library clean-objects clean-source clean-driver clean-makefile

clean: clean-library clean-objects clean-source
distclean: clean clean-driver clean-makefile
.PHONY: clean distclean
