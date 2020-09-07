//#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <RcppParallel.h>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <queue>

using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]
struct eucDist : public Worker {
   
   // input matrix to read from
   const RMatrix<double> mat;
      
   // output matrix to write to
   //RMatrix<double> rmat;
   RMatrix<double> output;

  // initialize from Rcpp input and output matrices (the RMatrix class
   // can be automatically converted to from the Rcpp matrix type)
   eucDist(const NumericMatrix mat, NumericMatrix output)
      : mat(mat), output(output) {}
   
   // function call operator that work for the specified range (begin/end)
   void operator()(std::size_t begin, std::size_t end) {
      
//      std::transform(mat.begin() + begin, 
//                     mat.begin() + end, 
//                     output.begin() + begin, 
//                     ::sqrt);

          
      for (std::size_t i = begin-1; i < end; i++) {
        for (std::size_t j = i+1; j < mat.nrow(); j++) {
            
            //rows we will operate on
            RMatrix<double>::Row row1 = mat.row(i);
            RMatrix<double>::Row row2 = mat.row(j);
            
            // calculate euclidean distance
           
        double rval = 0;
        for (unsigned int k = 0; k < row1.length(); ++k) {
          rval += (row1[k] - row2[k]) * (row1[k] - row2[k]);
        }
       
      rval=std::sqrt(rval);                      
      output(i,j) = rval;
      output(j,i) = rval;
         
      }
     }  
   }
};

// [[Rcpp::export]]
NumericMatrix eucDist_parallel(NumericMatrix mat) {
  
   // allocate the matrix we will return
   NumericMatrix output(mat.nrow(), mat.nrow());

   // create the worker
   eucDist eucdist(mat, output);
     
   // call it with parallelFor
   parallelFor(0, mat.nrow(), eucdist);

   return output;
}

