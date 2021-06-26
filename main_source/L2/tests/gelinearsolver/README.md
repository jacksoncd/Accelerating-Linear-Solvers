# Interface source code

## IpVitisSolverInterface.cpp/IpVitisSolverInterface.hpp 
Interface source and header files

## kernelgelinearsolver.cpp 
Kernel top level function. This calls the [solver](../../include/hw/LinearSolver/gelinearsolver.hpp)

## TutorialCpp_main.cpp/TutorialCpp_main.hpp 
Microbenchmark source and header files

## test_gelinearsolver.cpp  
Original test host provided by Xilinx

## IpAlgBuilder.cpp 
Modified IPOPT file so that IPOPT finds the interface when its linear solver option is set to custom. This is later changed on the IPOPT version installed on host machine.

## examples/
Folder of IPOPT example scripts

## .ini files
Files to be passed to the Vitis compiler for the specific devices

