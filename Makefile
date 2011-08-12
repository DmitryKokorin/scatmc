# Makefile
TARGET=$(shell basename `pwd`)
SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:%.cpp=%.o)


#debug
#CPPFLAGS	+=-D__DEBUG__ -g

#profiling
#CPPFLAGS	+= -pg
#LDFLAGS	+= -pg

#arch
CPPFLAGS	+=-march=native -O3 -ffast-math -pipe

#warnings
CPPFLAGS	+=-W -Wall -Weffc++ -Werror -pedantic

#experimental: Henyey-Greenstein phase function
#CPPFLAGS    +=-DEXPERIMENTAL=1

#openmp
CPPFLAGS	+=-fopenmp
LDFLAGS		+=-fopenmp


all: $(TARGET)

$(OBJECTS): $(SOURCES)

$(TARGET): $(OBJECTS)
	$(CXX) -o $(TARGET) $(LDFLAGS) $(OBJECTS) $(LOADLIBES) $(LDLIBS)
  
clean:
	$(RM) $(OBJECTS) $(TARGET)
  
.PHONY: all clean
