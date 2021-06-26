/*

Header file for the interface to the Vitis library solver
Base class is SparseSymLinearSolver

*/

#ifndef __IPVITISSOLVERINTERFACE_HPP_
#define __IPVITISSOLVERINTERFACE_HPP_

#include <iostream>
#include <string.h>
#include <sys/time.h>
#include <algorithm>
#include <cmath>
#include <fstream>

#include "xcl2.hpp"

#include "IpSparseSymLinearSolverInterface.hpp"

//#define MAXN 998

namespace Ipopt
{
class VitisSolverInterface: public SparseSymLinearSolverInterface
{
private:
  
  // Memory allocator
    template <typename T>
  T* aligned_alloc(std::size_t num) {
      void* ptr = nullptr;
      if (posix_memalign(&ptr, 4096, num * sizeof(T))) {
          throw std::bad_alloc();
      }
      return reinterpret_cast<T*>(ptr);
  }
  
  // Time difference
  unsigned long diff(const struct timeval* newTime, const struct timeval* oldTime) {
    return (newTime->tv_sec - oldTime->tv_sec) * 1000000 + (newTime->tv_usec - oldTime->tv_usec);
  }
  
  
  // Matrix variables
  int numneg_; // Number of negative eigenvalues
  double * val_; // Ptr for values
  
  Index matrix_nonzeros; // Number of non zeros values in A
  Index matrix_dimension; // Dimension of A
  Index num_rhs;     // The number of different RHS we are solving
  int dataA_size;    // Size of array A
  int dataB_size;    // Size of array B
  
  // OpenCL variables
  std::string xclbin_path; // Path for FPGA binary file
  cl::Device device; // Chosen device
  
  cl::Context context; // OpenCL context
  cl::CommandQueue q; // OpenCL q
  std::string devName; // Device name string
  
  cl::Program::Binaries xclBins; // OpenCL binaries
  cl::Program program; // OpenCL program
  cl::Kernel kernel_gelinearsolver_0; // Device kernel
  
  std::vector<cl::Device> devices; // Vector of devices
  

  
  // Time variables
  //struct timeval tstart, tinit_parse, tplatform_setup, tbuffer_setup, tbuffer_transfer1, tkernel_setup, tkernel_launch, tbuffer_transfer2; // Variables to measure time



public:


  // Constructor
  VitisSolverInterface():val_(nullptr),devices(1)
  {
  };
  
  static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
   );
  
  ~VitisSolverInterface();
  
  int SetBinaryPath(std::string binary_path);
  
  
  
  bool InitializeImpl(
      const OptionsList& options,
      const std::string& prefix
   );
   
   ESymSolverStatus InitializeStructure(
      Index        dim,
      Index        nonzeros,
      const Index* ia,
      const Index* ja
   );

   double* GetValuesArrayPtr();

   ESymSolverStatus MultiSolve(
      bool         new_matrix,
      const Index* ia,
      const Index* ja,
      Index        nrhs,
      double*      rhs_vals,
      bool         check_NegEVals,
      Index        numberOfNegEVals
   );
   
   

   // Retrive number of eigenvalues (unused)

   Index NumberOfNegEVals() const
   {
      return numneg_;
   }
   
   // IPOPT calling for increased quality of solution (unused)
   bool IncreaseQuality()
   {
     return false;
   }

   

   // If solver provides inertia

   bool ProvidesInertia() const
   {
      return false;
   }

   

   // Set Triplet format
   EMatrixFormat MatrixFormat() const
   {
      return Triplet_Format;
   }




   // If solver provides degeneracy detection
   bool ProvidesDegeneracyDetection() const
   {
      return false;
   }



   // If solver determines dependant rows
   ESymSolverStatus DetermineDependentRows(
      const Index*      /*ia*/,
      const Index*      /*ja*/,
      std::list<Index>& /*c_deps*/
   )
   {
      return SYMSOLVER_FATAL_ERROR;
   }
   
};

} //namespace ipopt


#endif
   
   
