#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

    
// [[Rcpp::export]]
NumericMatrix EvSync_sign_serial(const int i, NumericVector ESobs, NumericMatrix ESnull, NumericVector tlen, NumericVector vec, NumericVector cellIDs) {
  
   NumericMatrix output(ESnull.nrow()-1, 4);
   int lx = tlen[vec[i - 1] - 1];
   int count = 0;
     
     for (int j = i; j < ESobs.length(); ++j) {
    
          int ly = tlen[vec[0] + j - 1];
          int Qthresh;
          int Qobs;
          Qobs = ESobs[j];
          
             for (int k = 0; k < ESnull.nrow(); k++) {
                                                     
                  if(ESnull(k,0) == lx && ESnull(k,1) == ly) {
          
                    Qthresh = ESnull(k,2);
                    
                      if(Qobs > Qthresh) {
		                           
                         output(count, 0) = cellIDs[vec[i - 1] - 1]; 
			                   output(count, 1) = cellIDs[vec[0] + j - 1];
                         output(count, 2) = Qobs;                            
                         output(count, 3) = Qthresh;
                 
                         count++;
                                         
                      }                              

                 }
             }
        }


   return output;
}

