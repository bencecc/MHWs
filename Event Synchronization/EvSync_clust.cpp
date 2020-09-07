#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <RcppParallel.h>

using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

    
// [[Rcpp::export]]
int EvSync_compute(std::vector<int> ex, std::vector<int> ey) {
  
// compute event synchronization given two vectors of time indices of event occurrences
  
            int lx = ex.size();
            int ly = ey.size(); 
              
            int dst, tmp;
            double tau;
            int taumax = 10;
            int m, n; 
            
            int count;
            int Q;
            
            count = 0;
            
                for(m = 1; m < lx - 1; m++) {
	                
                  for(n = 1; n < ly - 1; n++) {
                                 
                      dst = ex[m] - ey[n];
		                  if(dst > taumax)
                            continue;
                      tmp = ex[m+1] - ex[m];
                      if(tmp > ex[m] - ex[m-1])
                          tmp = ex[m] - ex[m-1];
                      tau = ey[n+1] - ey[n];
                      if(tau > ey[n] - ey[n-1])
                          tau = ey[n] - ey[n-1];
                      if(tau > tmp)
	                        tau = tmp;
                      tau /= 2;
                      if(abs(ex[m] - ey[n]) <= taumax && abs(ex[m] - ey[n]) < tau)
                          count++;
                                  if(dst < -taumax)
						                                break;                                         
 
                   }
              }
               
         
           if(lx < 3 ) {
             
             Q = 0;
           
           }
           
           else if(ly < 3) {
             
             Q = 0;
           
           }
           
           else {
           
             Q = count;
           
           }   
          
          return Q;
  
}

struct EvSync : public Worker {
   
   // input matrix
   const RMatrix<double> mat;
   const RMatrix<double> refmat;
      
   // output matrix
   RMatrix<double> output;
   
   EvSync(const NumericMatrix mat, const NumericMatrix refmat, NumericMatrix output) : mat(mat), refmat(refmat), output(output) {}
   
   void operator()(std::size_t begin, std::size_t end) {
   
     for (std::size_t i = begin-1; i < end; i++) {
  
       for (std::size_t j = i+1; j < refmat.ncol(); j++) {

            RMatrix<double>::Column col1 = mat.column(i);    
            RMatrix<double>::Column col2 = refmat.column(j);  
                        
            std::vector<int>  ex(col1.length());
            std::vector<int>  ey(col2.length());
            int lx, ly;
            lx = 0;
            ly = 0;
            
            for (unsigned k = 0; k < col1.length(); ++k) {
                
                  if(col1[k] > 0) {
                  
                      ex[lx++] = k;
                      
                   }
                  
                  if (col2[k] > 0) {
                  
                      ey[ly++] = k;
                  
                  }
             }
 
               int obsQ;
               obsQ = EvSync_compute(ex,ey);
               output(j, i) = double(obsQ);
               
        }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix EvSync_clust(NumericMatrix mat, NumericMatrix refmat) {
  
   NumericMatrix output(refmat.ncol(), mat.ncol());

   // create the worker
   EvSync evsync(mat, refmat, output);
     
   // call it with parallelFor
   parallelFor(0, mat.ncol(), evsync);

   return output;
}

