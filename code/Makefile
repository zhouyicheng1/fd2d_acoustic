#******************************************************************************#
#*                 Makefile for FD2Dflatpsv-staggered package                 *#
#*                                                                            *#
#*  Author:                                                                   *#
#*      Wenzhong CAO    Email: caowz@mail.ustc.edu.cn                         *#
#*      Wei ZHANG       Email: zhangwei.zw@gmail.com                          *#
#*                                                                            *#
#*                                                                            *#
#*                 Copyright (C) 2017 Wenzhong CAO & Wei ZHANG                *#
#******************************************************************************#

# $Date: 2009-01-15 15:46:41 -0500 (Thu, 15 Jan 2009) $
# $Revision: 281 $
# $LastChangedBy: zhangw $

#######################################################################
#                 Compiler and macro flags setting                    #
#######################################################################

#WHEREAMI := w520
WHEREAMI := cwz
#WHEREAMI := pku
#WHEREAMI := ustc
#WHEREAMI := generic

# Please choose COMPILER in the included Makefile.opt.$(WHEREAMI)

#LITE := ON
#DEBUG := ON
#STATIC := ON
#GETARG := ON
#VERBOSE :=ON

#WITHQS := ON
#DataTypeDouble := ON

FourthORDER :=ON

DFLAG_LIST := LITE DEBUG STATIC GETARG VERBOSE \
              DataTypeDouble \
			  WITHQS FourthORDER 

DFLAGS := $(foreach flag,$(DFLAG_LIST),$(if $($(flag)),-D$(flag),)) $(DFLAGS)
DFLAGS := $(strip $(DFLAGS))

#######################################################################
#            directories, source files and target names               #
#######################################################################
skeldirs := OBJ bin src
OBJDIR := ./OBJ
BINDIR := ./bin
SRCDIR := ./src
FPPDIR := ./srcF
#FPPDIR := ./src.yxn

SRC_MOD  := mod_constants.f90 mod_math.f90 mod_string.f90 mod_nfseis.f90 \
		  mod_para.f90 mod_nompi.f90 mod_grid.f90 \
          mod_io.f90 mod_media.f90 mod_staggerfd.f90 mod_src.f90 mod_injec.f90 \
          mod_abs_exp.f90 mod_abs_npml_cfs.f90 mod_abs_npml_dpole.f90 
          #mod_interformetry.f90 mod_interformetry_solver.f90 mod_interformetry_io.f90
          
OBJ_MOD  :=  $(foreach file,$(SRC_MOD),$(OBJDIR)/$(file:.f90=.o))

SRC_GRID     := seis2d_grid.f90
EXE_GRID     := seis2d_grid
OBJ_GRID     := $(foreach file,$(SRC_GRID),$(OBJDIR)/$(file:.f90=.o))

SRC_MEDIA     := seis2d_media.f90
EXE_MEDIA     := seis2d_media
OBJ_MEDIA     := $(foreach file,$(SRC_MEDIA),$(OBJDIR)/$(file:.f90=.o))

SRC_SOURCE    := seis2d_source.f90
EXE_SOURCE    := seis2d_source
OBJ_SOURCE    := $(foreach file,$(SRC_SOURCE),$(OBJDIR)/$(file:.f90=.o))

SRC_STATION   := seis2d_station.f90
EXE_STATION   := seis2d_station
OBJ_STATION   := $(foreach file,$(SRC_STATION),$(OBJDIR)/$(file:.f90=.o))

SRC_WAVE_INJEC := seis2d_wave_injec.f90
EXE_WAVE_INJEC  := seis2d_wave_injec
OBJ_WAVE_INJEC  := $(foreach file,$(SRC_WAVE_INJEC),$(OBJDIR)/$(file:.f90=.o))

SRC_WAVE_PSV := seis2d_wave_psv.f90
EXE_WAVE_PSV := seis2d_wave_psv
OBJ_WAVE_PSV := $(foreach file,$(SRC_WAVE_PSV),$(OBJDIR)/$(file:.f90=.o))

#SRC_INT := seis2d_interformetry.f90
#EXE_INT := seis2d_interformetry
#OBJ_INT := $(foreach file,$(SRC_INT),$(OBJDIR)/$(file:.f90=.o))
vpath %.F90 $(FPPDIR)
vpath %.h   $(FPPDIR)

FPP := /usr/bin/cpp
FPPFLAGS := -P -traditional $(foreach dir,$(FPPDIR),-I$(dir)) $(DFLAGS) 

#######################################################################
#                     compiler and option                             #
#######################################################################
include Makefile.opt.$(WHEREAMI)

#######################################################################
#                            target                                   #
#######################################################################
phony_list := skel all preprocess solver postprocess util
.PHONY: $(phony_list)

all: skel preprocess solver 
#all: skel preprocess solver postprocess kernel plus

preprocess: $(BINDIR)/$(EXE_GRID)  \
            $(BINDIR)/$(EXE_MEDIA) \
            $(BINDIR)/$(EXE_SOURCE) $(BINDIR)/$(EXE_STATION)
solver:  $(BINDIR)/$(EXE_WAVE_INJEC) $(BINDIR)/$(EXE_WAVE_PSV)
kernel: $(BINDIR)/$(EXE_KERNEL) $(BINDIR)/$(EXE_KERNEL_MPI)
plus: $(BINDIR)/$(EXE_GATHER) $(BINDIR)/$(EXE_CONVERT) $(BINDIR)/$(EXE_GRTM) \
      $(BINDIR)/$(EXE_SNAP2SEIS)
skel:
	@mkdir -p $(skeldirs)

$(BINDIR)/$(EXE_WAVE_INJEC): $(OBJ_MOD)  $(OBJ_WAVE_INJEC)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_WAVE_INJEC) $(LDFLAGS)
$(BINDIR)/$(EXE_WAVE_PSV): $(OBJ_MOD)  $(OBJ_WAVE_PSV)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_WAVE_PSV) $(LDFLAGS)
$(BINDIR)/$(EXE_GRID): $(OBJ_MOD) $(OBJ_GRID)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_GRID) $(LDFLAGS)
$(BINDIR)/$(EXE_MEDIA): $(OBJ_MOD) $(OBJ_MEDIA)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_MEDIA) $(LDFLAGS)
$(BINDIR)/$(EXE_SOURCE): $(OBJ_MOD) $(OBJ_SOURCE)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_SOURCE) $(LDFLAGS)
$(BINDIR)/$(EXE_STATION): $(OBJ_MOD) $(OBJ_STATION)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_STATION) $(LDFLAGS)

RM := rm
cleanexe:
	$(RM) -f $(BINDIR)/*
cleanobj:
	$(RM) -f $(OBJDIR)/*
cleanf90:
	$(RM) -f $(SRCDIR)/*
cleanall: cleanexe cleanobj cleanf90
distclean: cleanexe cleanobj cleanf90

#######################################################################
#                        suffixes rules                               #
#######################################################################

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

%.f90 : %.F90
	$(FPP) $(FPPFLAGS) $< > $(SRCDIR)/$(@F)

$(OBJ_GRID) : $(SRC_GRID)
	$(FC) $(FFLAGSMP) $(SRCDIR)/$(<F) -o $@

$(OBJDIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(SRCDIR)/$(<F) -o $@

# vim:ft=make:ts=4:sw=4:nu:et:ai:
