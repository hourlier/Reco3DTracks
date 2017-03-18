
# NAME OF LIBRARY/PACKAGE
PACKAGE_NAME = Reco3DTracks

# ADD BINDARY NAMES HERE
PROGRAMS = run_Reco3D_larlite

PROGRAM_SOURCE = $(addsuffix .cxx, $(PROGRAMS))
SOURCES = $(filter-out $(PROGRAM_SOURCE), $(wildcard *.cxx))
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))

# Include your header file location
# include options for this package
INCFLAGS  = -g -I.                       #Include itself
INCFLAGS += $(shell larcv-config --includes) 
INCFLAGS += $(shell larcv-config --includes)/../app
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell larlite-config --includes)/../UserDev/BasicTool
INCFLAGS += $(shell larlite-config --includes)/../UserDev/SelectionTool
# example of adding in additional, non-core larlite packages
#INCFLAGS += $(shell larlite-config --includes)/../UserDev/BasicTool/GeoAlgo
#INCFLAGS += $(shell larlite-config --includes)/../UserDev/BasicTool/FhiclLite
INCFLAGS += $(shell larlitecv-config --includes)

CXXFLAGS += -I. $(shell root-config --cflags) -g

CXXFLAGS += $(shell larlite-config --includes)
CXXFLAGS += $(shell larlite-config --includes)/../UserDev
CXXFLAGS += $(shell larcv-config --includes)
CXXFLAGS += $(shell larcv-config --includes)/../app
CXXFLAGS += $(shell larlitecv-config --includes)
CXXFLAGS += $(shell larlitecv-config --includes)/../app

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

# Include your shared object lib location
LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlitecv-config --libs)
LDFLAGS += $(shell root-config --libs) -lPhysics -lMatrix -g

# call the common GNUmakefile
include $(LARLITECV_BASEDIR)/Makefile/GNUmakefile.CORE

pkg_build:
	@rm -rf $(PROGRAMS)
	@cp $(PROGRAMS).cxx $(PROGRAMS).cxx~
	@rm -rf tmp.txt
pkg_clean:

# Add your program below with a space after the previous one.
# This makefile compiles all binaries specified below.

all: $(PROGRAMS)

$(PROGRAMS):
	make --directory=../..
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cxx -o $@ $(CXXFLAGS) $(LDFLAGS)
	@rm -rf *.dSYM