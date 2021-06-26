CXX = g++
CC = gcc

IP_DIR = /home/jacksoncd/git/Ipopt/src
XFLIB_DIR = /home/jacksoncd/solver-acceleration/linear_solver_lib/solver

LIBFILE = libvitislibrary.so

CXXFLAGS += -I$(XILINX_XRT)/include -I$(XILINX_VIVADO)/include -std=c++14  -Wall -Wno-unknown-pragmas -Wno-unused-label -pthread
CXXLINKFLAGS += -L$(XILINX_XRT)/lib -lxilinxopencl -lpthread -lrt -Wno-unused-label -Wno-narrowing -DVERBOSE 
CXXLINKFLAGS += -L$(XILINX_VIVADO)/lnx64/tools/fpo_v7_0 -Wl,--as-needed -lgmp -lmpfr -lIp_floating_point_v7_0_bitacc_cmodel 

CXXFLAGS += -I$(XFLIB_DIR)/L2/tests/gelinearsolver
CXXFLAGS += -I$(XFLIB_DIR)/ext/xcl2
CXXFLAGS += -I$(IP_DIR)/Algorithm/LinearSolvers
CXXFLAGS += -I$(IP_DIR)/Common
CXXFLAGS += -I$(IP_DIR)/Algorithm
CXXFLAGS += -I$(IP_DIR)/Interfaces
CXXFLAGS += -I$(IP_DIR)/LinAlg

MAX_SIZE ?= 100

CXXFLAGS += -DMAXN=$(MAX_SIZE)

CXXFLAGS += -O2
CXXLINKFLAGS += -O2

.PHONY: shared
shared: $(LIBFILE)

.PHONY: clean
clean:
	rm $(LIBFILE) xcl2.o IpVitisSolverInterface.o


xcl2.o: xcl2.cpp
	$(CXX) -c xcl2.cpp -fPIC -o xcl2.o $(CXXFLAGS)

IpVitisSolverInterface.o: IpVitisSolverInterface.cpp
	$(CXX) -c IpVitisSolverInterface.cpp -fPIC -o IpVitisSolverInterface.o $(CXXFLAGS)

$(LIBFILE): IpVitisSolverInterface.o xcl2.o
	$(CXX) -shared -o $(LIBFILE) xcl2.o IpVitisSolverInterface.o $(CXXLINKFLAGS)
