// Copyright (C) 2009 International Business Machines.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Author:  Andreas Waechter               IBM    2009-04-02

// This file is part of the Ipopt tutorial.  It is a correct version
// of a C++ implementation of the coding exercise problem (in AMPL
// formulation):
//
// param n := 4;
//
// var x {1..n} <= 0, >= -1.5, := -0.5;
//
// minimize obj:
//   sum{i in 1..n} (x[i]-1)^2;
//
// subject to constr {i in 2..n-1}:
//   (x[i]^2+1.5*x[i]-i/n)*cos(x[i+1]) - x[i-1] = 0;
//
// The constant term "i/n" in the constraint is supposed to be input data
//

#include <cstdio>
#include <fstream>
#include <sys/time.h>
#include <chrono>
#include <string>

#include "IpIpoptApplication.hpp"
#include "TutorialCpp_nlp.hpp"

using namespace Ipopt;

// Compute time difference
unsigned long diff(const struct timeval* newTime, const struct timeval* oldTime) {
    return (newTime->tv_sec - oldTime->tv_sec) * 1000000 + (newTime->tv_usec - oldTime->tv_usec);
}

int main(
   int argv,
   char* argc[]
)
{
    /***********
    Store timing values
    **************/
    //struct timeval tstart, tend;
    
    int number_iterations = std::atoi(argc[3]);
    
    for(int iteration = 1; iteration <= number_iterations; iteration++)
    {   
    
    FILE* file = fopen(argc[1],"a");

   // Number of variables
   Index N = iteration*std::atoi(argc[4]);
   
   fprintf(file,"%d ",N);
   
   auto start = std::chrono::steady_clock::now();

   // constant terms in the constraints
   Number* a = new double[N - 2];
   for( Index i = 0; i < N - 2; i++ )
   {
      a[i] = (double(i + 2)) / (double) N;
   }

   // Create a new instance of your nlp
   //  (use a SmartPtr, not raw)
   SmartPtr<TNLP> mynlp = new TutorialCpp_NLP(N, a);

   // Create a new instance of IpoptApplication
   //  (use a SmartPtr, not raw)
   SmartPtr<IpoptApplication> app = new IpoptApplication();

   // Change some options
   // Note: The following choices are only examples, they might not be
   //       suitable for your optimization problem.
   app->Options()->SetNumericValue("tol", 1e-7);
   app->Options()->SetStringValue("mu_strategy", "adaptive");
   app->Options()->SetStringValue("linear_solver", argc[2]);
   //app->Options()->SetIntegerValue("max_iter", 10);

   // Intialize the IpoptApplication and process the options
   app->Initialize();

   // Ask Ipopt to solve the problem
   ApplicationReturnStatus status = app->OptimizeTNLP(mynlp);
      
   if( status == Solve_Succeeded )
   {
      printf("\n\n*** The problem solved!\n");
   }
   else
   {
      printf("\n\n*** The problem FAILED!\n");
   }

   // As the SmartPtrs go out of scope, the reference count
   // will be decremented and the objects will automatically
   // be deleted.

   // However, we created the Number array for a here and have to delete it
   delete[] a;

   auto end = std::chrono::steady_clock::now();
   
   fprintf(file, "%f \n", std::chrono::duration <double, std::milli> (end-start).count());
   
   fclose(file);
   
   } //for loop

   return 0;
}
