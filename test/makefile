# uncomment to have debug build
#DEBUG ?= 1

SRC_DIR := ../src
VPATH := $(SRC_DIR)

CC := mpic++
# If directly using the underlying compiler, e.g., g++, then switches included in the above
# wrapper compiler need to be explicitly added. For MPICH, this can be done by the commented
# out line below. For OpenMPI, replace -show with -showme:compile and/or -showme:link and
# remove the filter-out function
#CFLAGS := $(filter-out g++, $(shell mpic++ -show)) # for compiler other than g++, replace here

CFLAGS := -O3 -std=c++17 -Wall
ifeq ($(DEBUG),1)
  CFLAGS += -g3 -O0
endif

INCL := -I$(SRC_DIR)
SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp, %.o, $(notdir $(SRCS)))
EXEC := test

all: $(EXEC)
$(EXEC): $(OBJS)
	$(CC) -o $@ $(OBJS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

clean:	
	rm -f $(OBJS) $(EXEC) *.o
	
.PHONY : all clean
