#
# This is an example GNUmakefile for my packages
#

# specific names for this package
DICT  = Kazu_AdrienDict
SHLIB = libKazu_Adrien.so
SOURCES = $(filter-out $(DICT).cxx, $(wildcard *.cxx))
FMWK_HEADERS = LinkDef.h $(DICT).h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
OBJECTS = $(SOURCES:.cxx=.o)

# include options for this package
INCFLAGS  = -I.                       #Include itself
INCFLAGS += -I$(LARCV_BASEDIR)/../
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell larcv-config --includes)
INCFLAGS += $(shell larlitecv-config --includes)
INCFLAGS += $(shell larlitecv-config --includes)/../app
INCFLAGS += $(shell geo2d-config --includes)
INCFLAGS += $(shell root-config --cflags) -g

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

# call kernel specific compiler setup
include $(LARLITE_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlitecv-config --libs)
LDFLAGS += $(shell root-config --libs) -g

include $(LARLITE_BASEDIR)/Makefile/GNUmakefile.CORE


