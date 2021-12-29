	######################################################
	#	    Directory based Makefile		     #
	#      Makefile Created on Dec 04 2019 By 	     #
	#	  	Mohammad Reza Eidi     		     #
	######################################################
	#  Makefile by Make				     #
	#  Run file by make runexe			     #
	#  Debuge by ddd debuger by make debug		     #
	#  Clean objects and modules by make claen	     #
	#  Clean outputs and logs by make cleanout	     #
	#  Using MKL Libraries 				     #
	######################################################

# Detecting the OS
UNAME := $(shell uname)

# Compiler and Linker
FC	:= ifort

# Flags
FCFLAGS	:= -g -mcmodel=large
FCFLAGS	:= -debug -traceback -mcmodel=large
FCFLAGS	:= -O3 -qopt-report=5 -mcmodel=large
FCFLAGS	:= -O3 -qopenmp -qopt-report=5 -xhost
FCFLAGS	:= -O3 -qopenmp -xhost  
# FCFLAGS	:= -fast -qopenmp -xhost

# -qopenmp 
# -traceback
# -heap-arrays
# -check all
# -check bounds
# -check uninit
# -ftrapuv
# -debug all
# -gen-interfaces
# -warn interfaces
# -qopt-report=5
# -xhost


# (OMP_NUM_THREADS=8  and  export OMP_NUM_THREADS)  or (export OMP_NUM_THREADS=8)
# echo $OMP_NUM_THREADS

# Libraries
ifeq ($(UNAME), Linux)
LIBS    :=   -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl 
endif
ifeq ($(UNAME), Darwin)
LIBS	:=    -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_rt -lpthread -lm -ldl
endif

#The Target Binary Program
PROGNAME:= FFT_HHG

# The Directories
OUTDIR	:= out
PLOTDIR	:= plot
SRCDIR	:= src
COMPDIR	:= exe
BASHDIR	:= bash

# File extensions
SRCEXT	:= f90
OBJEXT	:= o
MODEXT	:= mod


PROG	:= $(COMPDIR)/$(PROGNAME)

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
SOURCES	:= $(sort $(wildcard $(SRCDIR)/*.$(SRCEXT)))
OBJECTS_UNSORTED := $(patsubst $(SRCDIR)/%,$(COMPDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))

#Reverse sort
OBJECTS := $(shell for w in $(OBJECTS_UNSORTED); do echo $$w; done | sort -r)

#Defauilt Make
all: resources $(PROG)

mall: all runall

#Remake
remake: call mall

#Copy Resources from Resources Directory to Target Directory
resources: directories
	@cp $(SRCDIR)/*.$(SRCEXT) $(COMPDIR)/
# with the use @ the following command is not printed

#Make the Directories
directories:
	@mkdir -p $(COMPDIR)

# run exe + do bash jobs
runall: runexe runbash

runexe: 
	./$(PROG)

runbash:
	./$(BASHDIR)/bash_harm.sh
	./$(BASHDIR)/bash_plot.sh
	
debug: 
	kdbg ./$(PROG)

# clean all 
call: clean cout cplot

# clean compiled file
clean: 
	@$(RM) -rf $(COMPDIR)
	
# clean outputs	
cout: 
	@$(RM) -f $(OUTDIR)/*
	
# clean plots
cplot: 
	@$(RM) -f $(PLOTDIR)/*

#Compile
$(COMPDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	$(FC) $(FCFLAGS) -o $@ -c $< 
	
$(PROG): $(OBJECTS) 
	$(FC) $(FCFLAGS) -o $@ $(OBJECTS) $(LIBS)
	@echo   built: `date "+%a %b %e %H:%M:%S %Z %Y"`> $(COMPDIR)/built_time.log
	@rm -f *.$(MODEXT)
	@rm -f $(COMPDIR)/*.$(OBJEXT)
	@rm -f $(COMPDIR)/*.$(SRCEXT)
	
# target: prerequisites
#	recipes (what to do)
# $@  is equal to the target
# $<  is equal to the first prerequisite : is not working if first prerequisite is an array
# $^  is equal to all prerequisites
	
# Non-File Targets
.PHONY: all remake runall runexe runbash runplot debug clean call mall cout cplot resources directories 
