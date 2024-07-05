###########################
# project directory paths #
###########################
PROJECT_ROOT := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))/..
rootdir = $(realpath .)

###############################################################
# Compiler specifications
# VPATH : Is the list of directories to be searched for missing
#     prerequisites. (Generally used for source files).
#     The best approach is to seperate the list by spaces.
###############################################################
VPATH = $(PROJECT_ROOT)/demo.2.mesh $(PROJECT_ROOT)/include

# Which compilers should be used
CXX = mpicxx
FLAGS = -fpermissive -fPIC -g

METIS_DIR = ${PROJECT_ROOT}/metis-5.0.3/install
METIS_INCLUDES = ${METIS_DIR}/include
METISLIB_PATH = -L${METIS_DIR}/lib -Wl,-rpath=${METIS_DIR}/lib
METISLIB = -lmetis
LIBS = -lm

INCLUDES = -I${METIS_INCLUDES}
# name of executable
EXE=MetisDemo

all: main.o 
	$(CXX) main.o $(LIBS) $(METISLIB_PATH) $(METISLIB) -o $(EXE)
	@echo "  >> main.cxx changed! Recompiled!\n"

main.o: main.cxx
	$(CXX) $(INCLUDES) -c main.cxx $(FLAGS)
	@echo "  >> main.cxx changed! Recompiled!\n"

clean:
	rm *.o *.vtu $(EXE)
	@echo "  >> Removed *.o files!\n"
