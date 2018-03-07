PDFINC ?= $(shell lhapdf-config --cflags)
APFELINC ?= $(shell apfel-config --cppflags)
GSLINC ?= $(shell gsl-config --cflags)


LDFLAGS = -L$(shell lhapdf-config --libdir) -lLHAPDF
LDFLAGS += -L$(shell apfel-config --libdir ) -lAPFEL
LDFLAGS += $(shell gsl-config --libs)


CXX    = c++
CXXCPP = c++

all: pdf

pdf: main.o
	$(CXX) main.o $(LDFLAGS) -o pdf 


main.o: main.C 
					$(CXXCPP) -fPIC $(PDFINC) $(APFELINC) $(GSLINC) -o main.o -c main.C



clean:
				rm -r pdf main.o 
