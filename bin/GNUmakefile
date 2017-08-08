# Include your header file location
CXXFLAGS += -I$(LARLITE_USERDEVDIR)/Kazu
CXXFLAGS += -I$(LARCV_BASEDIR)/../
CXXFLAGS += -I. $(shell root-config --cflags) -g
CXXFLAGS += $(shell larlite-config --includes)
CXXFLAGS += $(shell larcv-config --includes)
CXXFLAGS += $(shell larlitecv-config --includes)
CXXFLAGS += $(shell larlitecv-config --includes)/../app
CXXFLAGS += $(shell geo2d-config --includes)
CXXFLAGS += $(shell root-config --cflags) -g

# Include your shared object lib location
LDFLAGS += -L$(LARLITE_LIBDIR) -lKazu_Adrien 
LDFLAGS += $(shell root-config --libs) -lPhysics -lMatrix -g
LDFLAGS += $(shell larlite-config --libs)
LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlitecv-config --libs)
LDFLAGS += $(shell root-config --libs) -g

# platform-specific options
OSNAME = $(shell uname -s)
include $(LARLITE_BASEDIR)/Makefile/Makefile.${OSNAME}

# Add your program below with a space after the previous one.
# This makefile compiles all binaries specified below.
PROGRAMS = ReadReiFile StoppingPower

all:		$(PROGRAMS)

$(PROGRAMS):
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cc -o $@ $(CXXFLAGS) $(LDFLAGS)
	@rm -rf *.dSYM
clean:	
	rm -f $(PROGRAMS)
