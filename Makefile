# Define compilers
CXX := g++
CXXFLAGS := -Wall 
CXXPERFORMANCEFLAGS :=

# Define libraries and includes
LIBFLAGS := -L/usr/lib/lapack -llapack -L/usr/lib/libblas -lblas 
INCFLAGS := -I/usr/lib

# Other variables
EXTARGET := main

# Targets
.PHONY: all run clean $(EXTARGET)

all: $(EXTARGET)

OBJS = main.o

$(EXTARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LIBFLAGS) 

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $<

clean:
	-@rm -f $(EXTARGET) *.o *.a
