/*
 * Copyright 2019 Xilinx, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
 
#include <iostream>
#include <string.h>
#include <sys/time.h>
#include <algorithm>
#include <cmath>

#include "xcl2.hpp"

#define NUM_ROWS 500
#define NUM_ZEROS 7000
#define MAX_MATRIX_VALUE 500
#define NUM_ITERATIONS 1

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

    /**************
     Platform related operations
     *******************/
     
    printf("INFO: Running on HW \n");
      std::string xclbin_path = "/home/jacksoncd/solver-acceleration/linear_solver_lib/solver/L2/tests/gelinearsolver/build_dir.hw.xilinx_u50_gen3x16_xdma_201920_3/kernel_gelinearsolver.xclbin";
      
    std::vector<cl::Device> devices = xcl::get_xil_devices();
    cl::Device device = devices[0];

    // Creating Context and Command Queue for selected Device
    cl::Context context(device);
    cl::CommandQueue q(context, device, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);
    std::string devName = device.getInfo<CL_DEVICE_NAME>();
    printf("INFO: Found Device=%s\n", devName.c_str());

    cl::Program::Binaries xclBins = xcl::import_binary_file(xclbin_path);
    devices.resize(1);
    cl::Program program(context, devices, xclBins);
    cl::Kernel kernel_gelinearsolver_0(program, "kernel_gelinearsolver_0");
    std::cout << "INFO: Kernel has been created" << std::endl;
    
    // Create vector variables
    std::vector<std::vector<cl::Event>> kernel_evt(2);
    kernel_evt[0].resize(1);
    kernel_evt[1].resize(1);
    std::vector<cl::Buffer> buffer(2);


    /*********************
    Continuous Iterations
    ********************/
    
    // Timing variables
    struct timeval tstart, ttrans1, tlaunch, ttrans2;
    
    // Txt file open
    FILE* fp = fopen("matrix_values.txt","w");
    
    
     for(int iteration = 0; iteration <= 0; iteration++)
    {
    
    // Set debug mode and print to txt
    //int debug_mode = iteration;
    
    /*********************
    Sparse symmetric matrix intialisation
    *********************/
    
    // Matrix sizes and max value
    int num_rhs = 1;
    double divisor = 5;
    
    int num_rows = NUM_ROWS;
    int num_zeros = NUM_ZEROS;
    int max_value = MAX_MATRIX_VALUE;
    
    int matrix_size = num_rows*num_rows;
    int b_size = num_rows;
    
    
    // Number of nonzero values
    int num_nonzeros = matrix_size - num_zeros;
    
    // Reseed the rand function to ensure the same values each time
    srand(1);
        
    // Initialise matrix
    double * dataA;
    dataA = aligned_alloc<double>(matrix_size);
    double * dataB;
    dataB = aligned_alloc<double>(b_size);
    
    int ia[num_nonzeros];
    int ja[num_nonzeros];
    double matrix_values[num_nonzeros];
    
    for(int i = 0; i < num_nonzeros; i++)
    {
        ia[i] = rand() % num_rows;
        ja[i] = rand() % num_rows;
        matrix_values[i] = (rand() % max_value)/divisor;
    }
    
    
    // Populate the initial A array with zeros
     for(int i = 0; i < matrix_size; i++)
    {
       dataA[i] = 0;   
    }
    
    // Populate the array with randomly placed values
    for(int i = 0; i < num_nonzeros; i++)
        {   
            if(matrix_values[i] != 0)
            {
                dataA[num_rows*ia[i] + ja[i]] = matrix_values[i];
            }
     }
    
     
     // Ensure that matrix is symmetric by copying the lower triangle to the upper
     for(int i = 0; i < num_rows; i++)
     {
         for(int j = 0; j < num_rows; j++)
         {
             dataA[num_rows*i + j] = dataA[num_rows*j + i];
         }
     }
     
     // Initialise B
     for(int i = 0; i < b_size; i++)
     {
         dataB[i] = rand() % max_value;
     }
     
     // Print dataA and dataB to txt file
     for(int i = 0; i < num_rows; i++)
     {
         for(int j = 0; j < num_rows; j++)
         {
             fprintf(fp,"%f ", dataA[i*num_rows + j]);
         }
         
         fprintf(fp,"\n");
     }
     
     fprintf(fp,"\n \n");
     
     for(int i = 0; i < b_size; i++)
     {
         fprintf(fp, "%f ", dataB[i]);
     }
     
     fprintf(fp, "\n \n");
     
    

    gettimeofday(&tstart, 0);
    
    
    // Create device buffer and map dev buf to host buf
    buffer[0] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE,
                           sizeof(double) * matrix_size, dataA, NULL);
    buffer[1] = cl::Buffer(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_WRITE,
                           sizeof(double) * b_size, dataB, NULL);

    // Data transfer from host buffer to device buffer
    q.enqueueMigrateMemObjects({buffer[0],buffer[1]}, 0, nullptr, &kernel_evt[0][0]); // 0 : migrate from host to dev
    q.finish();
    
     gettimeofday(&ttrans1,0);
     
     int new_matrix = 1;
    

    // Setup kernel
    int debug_mode = 0;
    
    kernel_gelinearsolver_0.setArg(0, new_matrix);
    kernel_gelinearsolver_0.setArg(1, debug_mode);
    kernel_gelinearsolver_0.setArg(2, num_rhs);
    kernel_gelinearsolver_0.setArg(3, num_rows);
    kernel_gelinearsolver_0.setArg(4, buffer[0]);
    kernel_gelinearsolver_0.setArg(5, buffer[1]);
    q.finish();

    // Launch kernel and compute kernel execution time
    q.enqueueTask(kernel_gelinearsolver_0, nullptr, nullptr);
    q.finish();
    
     gettimeofday(&tlaunch,0);
    
    

    // Data transfer from device buffer to host buffer
    q.enqueueMigrateMemObjects({buffer[0], buffer[1]}, 1, nullptr, nullptr); // 1 : migrate from dev to host
    q.finish();
    
     gettimeofday(&ttrans2,0);

    for(int i = 0; i < b_size; i++)
    {
        if(std::isnan(dataB[i]))
        {
            printf("INFO : Matrix singular \n");
            break;
        }
    } 
    
    for(int i = 0; i < b_size; i++)
     {
         fprintf(fp, "%f ", dataB[i]);
     }
     
    
    
    
    
    
    free(dataA);
    free(dataB);
    
    // Output the time to txt file
    int trans1 = diff(&ttrans1,&tstart);
    int launch = diff(&tlaunch,&ttrans1);
    int trans2 = diff(&ttrans2,&tlaunch);

    //fprintf(fp,"Matrix dimension : %d \n",num_rows);
    //fprintf(fp,"First transfer : %d \n", trans1);
    //fprintf(fp,"Launch : %d \n", launch);
    //fprintf(fp,"Second transfer : %d \n", trans2);
    
    
    } // for loop
    
    
    // Txt file close
    fclose(fp);
    
    return 0;
}
