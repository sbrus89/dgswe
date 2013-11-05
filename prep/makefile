SHELL:=/bin/sh
#
#  Makefile to Build PADCIRC and its pre/post-processor 
#  University of Texas's parallel version of the Hydrodynamics
#  Modeling program ADCIRC by J.J. Westerink and R.A. Luettich
#  vjp 3/29/2001
#  updated JJW 3/17/2003
#  updated for 3D JGF spring 2003

########################################################################
#  Get Canonical Machine NAME from config.guess
#
NAME     := $(shell ./config.guess)
LIST     := $(subst -, ,$(NAME))
MACHINE  := $(word 1, $(LIST))
VENDOR   := $(word 2, $(LIST))
OS       := $(subst  $(MACHINE)-$(VENDOR)-,,$(strip $(NAME)))

include cmplrflags.mk

######################## Target Specific Rules ###################################

#                                      adcprep   
ifeq ($(BUILDTYPE),adcprep)
  CF:= $(PPFC)
  O_DIR:=odir1/
  FFLAGS:= $(FFLAGS1) $(DPRE) $(IMODS)$(O_DIR) 
  VPATH :=  ./ 
  PMSG_OBJ:=
endif

$(O_DIR):
	mkdir -p $@

######################### Module Source, Object, & Mod Files ######################

PREP_MSRC  =  presizes.F 
PGLO_MSRC  =  pre_global.F  

PREP_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PREP_MSRC) )
PGLO_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PGLO_MSRC) )

############################# Source & Object Files ##############################

METIS1_SRC =  metis.F
METIS2_SRC =  metis2.F
PREP_SRC   =  adcprep.F  decomp.F   read_global.F   prep.F   interp.F  machdep.F

METIS1_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(METIS1_SRC) )
METIS2_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(METIS2_SRC) )
PREP_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PREP_SRC) )

######################## compilation rules #####################################

$(O_DIR)%.o  : %.F
#  @echo depend $<
#  @echo target $@
	$(CF) -c $(FFLAGS) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi

#Casey 121126: Changed the f90 files to F90 extension to avoid conflict with SWAN.
$(O_DIR)%.o  : %.F90
	$(CF) -c $(FFLAGS) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi

########################## Executable Targets ##################################

.PHONY: all metis adcprep 

all :  metis adcprep 

ifeq ($(MAKELEVEL),0)
   metis:
	$(MAKE) -C ../metis/Lib/ CC="$(CC)"  CFLAGS="$(CFLAGS)"
   graphchk:
	$(MAKE) -C ../metis/Programs/ CC="$(CC)"  CFLAGS="$(CFLAGS)" 
	mv ../metis/graphchk ../work
   adcprep:
	$(MAKE) BUILDTYPE=adcprep  $@            
   adcprep2:
	$(MAKE) BUILDTYPE=adcprep2 $@ 
else
   adcprep::  $(O_DIR)
   adcprep2:: $(O_DIR)
   
   adcprep ::  $(METIS1_OBJ) $(PREP_OBJ) 
	$(CF) $(FFLAGS) -o $@  $(O_DIR)*.o  $(LIBS) $(LDFLAGS)

   adcprep2 :: $(PMSG_OBJ) $(METIS2_OBJ) $(PREP_OBJ) 
	$(CF) $(FFLAGS) -o $@  $(O_DIR)*.o  $(LIBS) $(MSGLIBS) $(LDFLAGS)

endif

########################## Misc Commands ####################################

clean:
	 rm -f odir*/*.o  odir*/*.mod 
clobber:
	rm -r -f odir*
	rm -f  graphchk adcprep adcprep2 adcpost adcirc padcirc adcswan padcswan \
                ../metis/Lib/*.o  ../metis/libmetis.a ../metis/Programs/*.o
help:
	@echo "This makefile supports the following:"
	@echo "make all      - makes all six targets"
	@echo "make adcprep  - makes the adcprep  executable"
	@echo "make adcprep2 - makes the adcprep2  executable"


########################## Defining the DAG  ####################################

#  adcprep & adcprep2 modules

$(O_DIR)presize.o     :  presize.F
$(O_DIR)pre_global.o  :  pre_global.F  $(PREP_MOBJ)

#  adcprep & adcprep2 

$(O_DIR)adcprep.o     :  adcprep.F  $(PGLO_MOBJ) $(PMSG_OBJ)
$(O_DIR)decomp.o      :  decomp.F   $(PGLO_MOBJ)
$(O_DIR)read_global.o :  read_global.F  $(PGLO_MOBJ)
$(O_DIR)prep.o        :  prep.F   $(PGLO_MOBJ)
$(O_DIR)interp.o      :  interp.F  
$(O_DIR)machdep.o     :  machdep.F
$(O_DIR)metis.o       :  metis.F $(PGLO_MOBJ)
$(O_DIR)metis2.o      :  metis2.F $(PGLO_MOBJ)
$(O_DIR)picomsg.o     :  picomsg.F

# graphchk

$(O_DIR)io.o		:  ../metis/Lib
