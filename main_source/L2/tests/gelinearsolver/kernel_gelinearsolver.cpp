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

#include "xf_solver_L2.hpp"
//#define NCU 36
//#define MAXN 998


extern "C" void kernel_gelinearsolver_0(int num_nonzeros, int new_matrix, int n, int num_rhs, int* A_rows, int* A_cols, double* A_vals, double* dataB) {
#pragma HLS INTERFACE m_axi port = A_rows bundle = gmem0 offset = slave max_read_burst_length = 128
#pragma HLS INTERFACE m_axi port = A_cols bundle = gmem0 offset = slave max_read_burst_length = 128
#pragma HLS INTERFACE m_axi port = A_vals bundle = gmem0 offset = slave max_read_burst_length = 128
#pragma HLS INTERFACE m_axi port = dataB bundle = gmem1 offset = slave max_read_burst_length = 128 max_write_burst_length = 128

#pragma HLS INTERFACE s_axilite port = num_nonzeros bundle = control
#pragma HLS INTERFACE s_axilite port = new_matrix bundle = control
#pragma HLS INTERFACE s_axilite port = n bundle = control
#pragma HLS INTERFACE s_axilite port = num_rhs bundle = control
#pragma HLS INTERFACE s_axilite port = A_rows bundle = control
#pragma HLS INTERFACE s_axilite port = A_cols bundle = control
#pragma HLS INTERFACE s_axilite port = A_vals bundle = control
#pragma HLS INTERFACE s_axilite port = dataB bundle = control
#pragma HLS INTERFACE s_axilite port = return bundle = control

    
    //int debug_mode = 0;
  
       // General linear solver
       xf::solver::gelinearsolver<double, MAXN_IN, NCU_IN>(num_nonzeros, new_matrix, n, num_rhs, A_rows, A_cols, A_vals, dataB);
    
}
