# Accelerating IPOPT linear solver using an FPGA
IPOPT is an open source nonlinear optimiser which requires a linear solver as part of its solving process. This repo contains the source code for a framework running this linear solver on an FPGA platform - specifically a Xilinx Alveo U50 card. The linear solver source code is modified from the Xilinx open source [Vitis Solver Library](https://github.com/Xilinx/Vitis_Libraries/tree/master/solver) to be optimised for use with IPOPT.


## Main Project
### [Main Source Code](main_source)
The main source code for the project. Both the IPOPT interface and main linear solver kernel source code:
#### [Kernel Code](main_source/L2/include/hw)
The C++ kernel source code split across two cpp files. [gelinearsolver.cpp](main_source/L2/include/hw/LinearSolver/gelinearsolver.hpp) is called by the top level [kernel](main_source/L2/tests/gelinearsovler/kernel_gelinearsolver.cpp) function and contains the solving algorithm. [getrf.hpp](main_source/L2/include/hw/MatrixDecomposition/getrf.hpp) contains the decomposition algorithm.

#### [Main source directory](main_source/L2/tests/gelinearsolver)
Main directory containing the top level kernel function for compilation. The interface source code and other platform related ini files are located here.

### [IPOPT Vitis Interface Library File](library_file_for_IPOPT)
The compile folder which creates a library file from the interface source code (symbolically linked to this directory from the [main source directory](main_source/L2/tests/gelinearsolver)). This is used by IPOPT when running with the Vitis Library Solver.

## Dependancies
### IPOPT
Instructions to install IPOPT are available at https://coin-or.github.io/Ipopt/. The source code is available at https://github.com/coin-or/Ipopt. This project used version 3.13.3.

Prior to running the configure script and install, the `IpAlgBuilder.cpp` file must be changed to include the Vitis interface. A copy of the updated cpp file is available in the main directory [here](main_source/L2/tests/gelinearsovler/IpAlgBuilder.cpp). [`IpVitisSolverInterface.hpp`](main_source/L2/tests/gelinearsovler/IpVitisSolverInterface.hpp) must also be symbolically linked into the `src/Algorithm/LinearSolvers` directory within IPOPT.

Usually IPOPT requires a linear solver package to be installed prior to its install and will therefore search for HSL MA57 by default when running the configuration script. Therefore to build IPOPT, it maybe necessary to run the configuration script with the `--without-hsl` flag should the compile fail due to a missing linear solver package.

### Vitis Software Package and XRT
Using the Alveo U50 FPGA requires the XRT library and the Vitis Software Platform. Both can be installed from the Xilinx website [here](https://www.xilinx.com/products/boards-and-kits/alveo/u50.html#gettingStarted). This project used version 2020.1.


## Using the code with IPOPT and Alveo U50 FPGA
### Xclbin Compile
To compile the xclbin file go to the [main source directory](main_source/L2/tests/gelinearsovler) and run:

`make xclbin`

This will compile the xclbin file and store it in the directory `build_dir.hw.xilinx_u50_gen3x16_xdma_201920_3`.

### Library file compile
IPOPT requires a library file, stored in `\usr\local\lib` which it calls when the Vits Library Solver is selected. First the location of the xclbin file must be set in the interface at line 73 of the [interface cpp file](main_source/L2/tests/gelinearsolver/IpVitisSolverInterface.cpp). This depends on the directory structure.
To compile go to the [library compile directory](library_file_for_IPOPT) and run:

`make shared`

This will produce a `.so` file which should be symbolically linked to `\usr\local\lib` for use with IPOPT. IPOPT has an environment variable that must be set to find the created library file, this can be set as follows:

`export LIBS=-lvitislibrary`

### Use
Set the linear solver IPOPT uses to custom. This can be done in an `ipopt.opt` file as follows:

`linear_solver custom`

IPOPT can then be run as usual and it will use the FPGA based solver.
