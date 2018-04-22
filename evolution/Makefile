PDFLIB   ?= $(shell lhapdf-config --libs)
PDFINC   ?= -I$(shell lhapdf-config --incdir)

APFELINC ?= $(shell apfel-config --cxxflags)
APFELLIB ?= $(shell apfel-config --ldflags)

GSLLIB   ?= $(shell gsl-config --libs)
GSLINC   ?= $(shell gsl-config --cflags)

CXX    = c++

CXXFLAGS =
LDFLAGS  =

all: test

test: test.o
	$(CXX) test.o $(APFELLIB) $(PDFLIB) $(GSLLIB) $(LDFLAGS) -o test

test.o: main.cpp
	$(CXX) -fPIC -std=c++14 $(PDFINC) $(APFELINC) $(GSLINC) -o test.o -c main.cpp

print-%  : ; @echo $* = $($*)

clean:
	rm -r test test.o
