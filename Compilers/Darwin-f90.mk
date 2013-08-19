# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2013 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for Absoft compiler on Mac (Darwin)
# -------------------------------------------------------------------------
#
# ARPACK_LIBDIR  ARPACK libary directory
# FC             Name of the fortran compiler to use
# FFLAGS         Flags to the fortran compiler
# CPP            Name of the C-preprocessor
# CPPFLAGS       Flags to the C-preprocessor
# CLEAN          Name of cleaning executable after C-preprocessing
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF libary directory
# LD             Program to load the objects into an executable
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#
               FC := f90
           FFLAGS := -N11 -YEXT_NAMES=LCS -YEXT_SFX=_ -YCFRL=1 -w
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional
          LDFLAGS :=
               AR := ar
          ARFLAGS := r
            MKDIR := mkdir -p
               RM := rm -f
           RANLIB := ranlib
             PERL := perl
             TEST := test

        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(SCRATCH_DIR)

#
# Library locations, can be overridden by environment variables.
#

ifdef USE_NETCDF4
        NC_CONFIG ?= nc-config
    NETCDF_INCDIR ?= $(shell $(NC_CONFIG) --prefix)/include
             LIBS := $(shell $(NC_CONFIG) --flibs)
else
    NETCDF_INCDIR ?= /usr/local/include
    NETCDF_LIBDIR ?= /usr/local/lib
             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf
endif

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef USE_MPI
         CPPFLAGS += -DMPI
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
endif

ifdef USE_DEBUG
           FFLAGS += -g
else
           FFLAGS += -O3
endif

ifdef USE_MCT
       MCT_INCDIR ?= /usr/local/mct/include
       MCT_LIBDIR ?= /usr/local/mct/lib
           FFLAGS += -I$(MCT_INCDIR)
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
endif

ifdef USE_ESMF
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
                     include $(ESMF_MK_DIR)/esmf.mk
           FFLAGS += $(ESMF_F90COMPILEPATHS)
             LIBS += $(ESMF_F90LINKPATHS) -lesmf -lC
endif

#
# Use full path of compiler.
#
               FC := $(shell which ${FC})
               LD := $(FC)

#
# Set free form format in source files to allow long string for
# local directory and compilation flags inside the code.
#

$(SCRATCH_DIR)/mod_ncparam.o: FFLAGS += -ffree
$(SCRATCH_DIR)/mod_strings.o: FFLAGS += -ffree
$(SCRATCH_DIR)/analytical.o: FFLAGS += -ffree
$(SCRATCH_DIR)/biology.o: FFLAGS += -ffree
ifdef USE_ADJOINT
$(SCRATCH_DIR)/ad_biology.o: FFLAGS += -ffree
endif
ifdef USE_REPRESENTER
$(SCRATCH_DIR)/rp_biology.o: FFLAGS += -ffree
endif
ifdef USE_TANGENT
$(SCRATCH_DIR)/tl_biology.o: FFLAGS += -ffree
endif

#
# Supress free format in SWAN source files since there are comments
# beyond column 72.
#

ifdef USE_SWAN

$(SCRATCH_DIR)/ocpcre.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/ocpids.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/ocpmix.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swancom1.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swancom2.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swancom3.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swancom4.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swancom5.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swanmain.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swanout1.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swanout2.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swanparll.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swanpre1.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swanpre2.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swanser.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swmod1.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/swmod2.o: FFLAGS += -ffixed
$(SCRATCH_DIR)/m_constants.o: FFLAGS += -ffree
$(SCRATCH_DIR)/m_fileio.o: FFLAGS += -ffree
$(SCRATCH_DIR)/mod_xnl4v5.o: FFLAGS += -ffree
$(SCRATCH_DIR)/serv_xnl4v5.o: FFLAGS += -ffree

endif
