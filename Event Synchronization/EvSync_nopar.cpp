#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int EvSync_nopar(NumericVector ex, NumericVector ey) {
  
  
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
