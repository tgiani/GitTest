PDFINC ?= $(shell lhapdf-config --cflags)
APFELINC ?= $(shell apfel-config --cppflags)
GSLINC ?= $(shell gsl-config --cflags)


LDFLAGS = -L$(shell lhapdf-config --libdir) -lLHAPDF
LDFLAGS += -L$(shell apfel-config --libdir ) -lAPFEL
LDFLAGS += $(shell gsl-config --libs)
LDFLAGS += -L. -lbbh
LDFLAGS += -L/usr/local/gfortran/lib -lgfortran



CXX    = c++
CXXCPP = c++

all: pdf

pdf: helpers.o main.o
	$(CXX) helpers.o main.o $(LDFLAGS) -o pdf 


helpers.o: helpers.C helpers.h
					$(CXXCPP) -fPIC $(PDFINC) $(APFELINC) $(GSLINC) -o helpers.o -c helpers.C

main.o: main.C 
					$(CXXCPP) -fPIC $(PDFINC) $(APFELINC) $(GSLINC) -o main.o -c main.C



clean:
				rm -r pdf main.o helpers.o
