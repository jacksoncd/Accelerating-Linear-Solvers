/*

Testing the interface

*/

#include "IpVitisSolverInterface.hpp"
using namespace Ipopt;

#include <fstream>


// Memory alignment
template <typename T>
T* aligned_alloc(std::size_t num) {
    void* ptr = nullptr;
    if (posix_memalign(&ptr, 4096, num * sizeof(T))) {
        throw std::bad_alloc();
    }
    return reinterpret_cast<T*>(ptr);
}

// Compute time difference
unsigned long diff(const struct timeval* newTime, const struct timeval* oldTime) {
    return (newTime->tv_sec - oldTime->tv_sec) * 1000000 + (newTime->tv_usec - oldTime->tv_usec);
}

// Arguments parser
class ArgParser {
   public:
    ArgParser(int& argc, const char** argv) {
        for (int i = 1; i < argc; ++i) mTokens.push_back(std::string(argv[i]));
    }
    bool getCmdOption(const std::string option, std::string& value) const {
        std::vector<std::string>::const_iterator itr;
        itr = std::find(this->mTokens.begin(), this->mTokens.end(), option);
        if (itr != this->mTokens.end() && ++itr != this->mTokens.end()) {
            value = *itr;
            return true;
        }
        return false;
    }

   private:
    std::vector<std::string> mTokens;
};

//! Core function of Cholesky benchmark
int main(int argc, const char* argv[]) {

   
    /********************
    Test IPOPT interface
    *********************/
    
    // Create instance
    VitisSolverInterface * solver_interface = new VitisSolverInterface;
    
    // Call initialiser
    const OptionsList options;
    const std::string prefix;
    solver_interface->InitializeImpl(options,prefix);
    
    // Initialise the structure using data from the HSL example
    Index dimension = 5;
    Index non_zeros = 8;
    
    // Data A values
    const Index ja[non_zeros] = {1,2,2,3,5,3,4,5};
    const Index * ja_ptr = ja;
    
    const Index ia[non_zeros] = {1,1,2,2,2,3,3,5};
    const Index * ia_ptr = ia;
    
    double nonzero_values[non_zeros] = {2.,1.,4.,1.,1.,3.,2.,2.};
    
    /***********
    Run and with different parts of kernel running to determine bottleneck
    ***************/
      
    
      solver_interface->InitializeStructure(dimension,non_zeros,ia_ptr,ja_ptr);
      
      
      // Find the location of the values array and populate
      double* value_pointer = solver_interface->GetValuesArrayPtr();
      
      for(int i = 0; i < non_zeros; i++){
          value_pointer[i] = nonzero_values[i];
      }
      
      // Allocate the values of B
      Index nrhs = 1;
      Index rhs_values_size = dimension*nrhs;
      double rhs_values[rhs_values_size] = {4.,12.,10.,4.,4.};
      double * rhs_values_ptr = new double[rhs_values_size];
      
      for(int i = 0; i < rhs_values_size; i++)
      {
          rhs_values_ptr[i] = rhs_values[i]; 
      }
      
      
      // Call the solve class method
      bool new_matrix = true;
      bool check_eigenvalues = false;
      Index eigenvalues = 0;
      
      ESymSolverStatus result = solver_interface->MultiSolve(new_matrix,ia_ptr,ja_ptr,nrhs,rhs_values_ptr,check_eigenvalues,eigenvalues);
      
      printf("Result : %d \n", result);
      
      // Print the result
      for(int i = 0; i < rhs_values_size; i++)
      {
          printf("x %d : %f \n",i,rhs_values_ptr[i]);
          //myfile << rhs_values_ptr[i] << " ";
      }
      
      //myfile << "\n";
      
      delete[] rhs_values_ptr;
    
   
   delete solver_interface;
    
}