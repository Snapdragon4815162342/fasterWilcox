#include <Rcpp.h>
#include <vector>
#include "sparse.h"

#include <iostream>
#include <algorithm>
#include <iterator>



// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>

#include <numeric>

#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame sparse_fasterWilcox_new(Rcpp::S4 A, NumericVector rowIndex, NumericVector col_id, NumericVector ref_col_id, bool verbose = 0, bool paired = 0, std::string alternative = "two.sided") {
  
  MySparseMatrix M =  MySparseMatrix(A);
  //R indices starts from 1
  // for (size_t i=0; i<rowIndex.size();i++)
  //   rowIndex[i]--;
  // for (size_t i=0; i<col_id.size();i++)
  //   col_id[i]--;
  // for (size_t i=0; i<ref_col_id.size();i++)
  //   ref_col_id[i]--;  
  
  // std::cout<<"\nM.x: ";
  // for(size_t i=0; i<M.x.size(); i++)
  //   std::cout<<M.x[i]<<" ";
  // 
  // std::cout<<"\nM.p: ";
  // for(size_t i=0; i<M.p.size(); i++)
  //   std::cout<<M.p[i]<<" ";
  // 
  // std::cout<<"\nM.i: ";
  // for(size_t i=0; i<M.i.size(); i++)
  //   std::cout<<M.i[i]<<" ";
  
  
  Rcpp::DataFrame results = M.sparseWilcox(rowIndex,
                 col_id,
                 ref_col_id,
                 paired,
                 alternative,
                 verbose);
  
  return results;
  
  }


