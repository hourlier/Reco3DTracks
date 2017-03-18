CXXFLAGS += $(shell root-config --cflags) -g
CXXFLAGS += $(shell larlite-config --includes)
CXXFLAGS += $(shell larlite-config --includes)/../UserDev
CXXFLAGS += $(shell larlite-config --includes)/../UserDev/BasicTool
CXXFLAGS += $(shell larlite-config --includes)/../UserDev/SelectionTool
CXXFLAGS += $(shell larcv-config --includes)
CXXFLAGS += $(shell larcv-config --includes)/../app
CXXFLAGS += $(shell larlitecv-config --includes)
#CXXFLAGS += $(shell larlitecv-config --includes)/../app

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

# Include your shared object lib location
LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlitecv-config --libs) -llarlitecv
LDFLAGS += $(shell root-config --libs) -lPhysics -lMatrix -g

# platform-specific options
OSNAME = $(shell uname -s)
include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

# Add your program below with a space after the previous one.
# This makefile compiles all binaries specified below.
PROGRAMS = run_Reco3D run_Reco3D_larlite

all:		$(PROGRAMS)

$(PROGRAMS):
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cc -o $@ $(CXXFLAGS) $(LDFLAGS)
	@rm -rf *.dSYM
clean:	
	rm -f $(PROGRAMS)
