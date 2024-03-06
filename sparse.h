#include <Rcpp.h>
#include <RcppCommon.h>
//[[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(cpp11)]]

// namespace MySparse {
//   class MySparseMatrix;
// }

// namespace MySparse {
class MySparseMatrix {
public:
  // fields
  Rcpp::NumericVector x;
  Rcpp::IntegerVector i, p, Dim;

  // constructor
  MySparseMatrix(const Rcpp::S4& s) {
    if (!s.hasSlot("x") || !s.hasSlot("p") || !s.hasSlot("i") || !s.hasSlot("Dim"))
      throw std::invalid_argument("Invalid S4 object, it must be a dgCMatrix");
    x = s.slot("x");
    i = s.slot("i");
    p = s.slot("p");
    Dim = s.slot("Dim");
  }

  int nrow() { return Dim[0]; };
  int ncol() { return Dim[1]; };

  Rcpp::NumericVector rowSumsColIndex(Rcpp::NumericVector colIndex) {
    Rcpp::NumericVector sums(Dim[0]);
    const int NCOL = colIndex.size();
    int col;
    for (int jcol = 0; jcol < NCOL; ++jcol){
      col = colIndex[jcol] - 1; // remove 1 since R indices (input) are 1-based
      for (int j = p[col]; j < p[col + 1]; ++j){
        sums[i[j]] += x[j];
      }
    }
    return sums;
  } // method rowSumsColIndex


  Rcpp::NumericVector sparseCountElementByRow_ColInd(Rcpp::NumericVector colIndex, double count_thr) {
    Rcpp::NumericVector sums(Dim[0]);
    const int NCOL = colIndex.size();
    int col;
    for (int jcol = 0; jcol < NCOL; ++jcol){
      col = colIndex[jcol] - 1; // remove 1 since R indices (input) are 1-based
      for (int j = p[col]; j < p[col + 1]; ++j){
        if(x[j] > count_thr){
          sums[i[j]] += 1;
        }

      }
    }
    return sums;
  } // method rowSumsColIndex


  Rcpp::NumericVector sparseCountElementByRow_RowInd_ColInd(Rcpp::NumericVector rowIndex, Rcpp::NumericVector colIndex, double count_thr) {

    Rcpp::NumericVector sums(Dim[0]);
    const int NCOL = colIndex.size();
    int col;
    for (int jcol = 0; jcol < NCOL; ++jcol){
      col = colIndex[jcol] - 1; // remove 1 since R indices (input) are 1-based
      for (int j = p[col]; j < p[col + 1]; ++j){
        if(x[j] > count_thr){
          sums[i[j]] += 1;
        }

      }
    }
    return sums[rowIndex -1];
  } // method sparseCountElementByRow_RowInd_ColInd



  /////////////////////////

  // DataFrame sparseWilc(Rcpp::NumericVector rowIndex,
  //                                Rcpp::NumericVector col_id,
  //                                Rcpp::NumericVector ref_col_id) {
  //
  //   Rcpp::NumericVector pvalues(pMat->nrow(), 1.0);
  //   Rcpp::NumericVector sign(pMat->nrow());
  //   sign.fill(NA_REAL); // NA_REAL is a constant representing NaN in Rcpp
  //
  //
  // }


  //////////////////////////

  Rcpp::NumericVector rowMeansRowColIndex(Rcpp::NumericVector rowIndex, Rcpp::NumericVector colIndex) {
    Rcpp::NumericVector means = rowMeansColIndex(colIndex);
    Rcpp::NumericVector row_means = means[rowIndex-1];
    return row_means;
  } // method rowMeansRowColIndex


  double MeanSparse() {
    double Mean = std::accumulate(x.begin(), x.end(), 0.0);
    return Mean/(Dim[0]*Dim[1]);
  } // method MeanSparse


  Rcpp::NumericVector rowMeansColIndex(Rcpp::NumericVector colIndex) {
    Rcpp::NumericVector sums = rowSumsColIndex(colIndex) / colIndex.size();
    return sums;
  } // method rowMeansColIndex


  double MeanSparseColIndex(Rcpp::NumericVector colIndex) {
    Rcpp::NumericVector geneMeans = rowMeansColIndex(colIndex);
    double mean = ( std::accumulate(geneMeans.begin(), geneMeans.end(), 0.0) ) / geneMeans.length();
    return mean;
  } // method MeanSparseColIndex

  Rcpp::List BigBootstrap(Rcpp::NumericVector colIndex,
                          Rcpp::NumericVector H0_genes,
                          Rcpp::NumericVector ligands,
                          Rcpp::NumericVector receptors){

    const int N_COL = colIndex.length();
    const int N_GENES = H0_genes.length();
    const int N_LIG = ligands.length();
    const int N_REC = receptors.length();

    // compute the mean of all genes
    Rcpp::NumericVector geneMeans = rowMeansColIndex(colIndex);

    // compute ligands and receptors means (subset of all genes mean)
    Rcpp::NumericVector ligands_avg = geneMeans[ligands-1];
    Rcpp::NumericVector receptors_avg = geneMeans[receptors-1];

    // compute H0 genes means (subset of all genes mean) to be used in matrix mean computation
    Rcpp::NumericVector H0geneMeans = geneMeans[H0_genes-1];
    // compute matrix mean as the mean of genes means
    double H0_mean = ( std::accumulate(H0geneMeans.begin(), H0geneMeans.end(), 0.0) ) / N_GENES;


    double H0_sd = 0.0;
    double K = H0_mean; // shift parameter
    double n = N_GENES * N_COL; // total number of elements
    Rcpp::NumericVector exv(Dim[0]);
    Rcpp::NumericVector ex2v(Dim[0]);
    Rcpp::NumericVector count_zero_by_row(Dim[0], N_COL*1.0);

    int col;
    for (int jcol = 0; jcol < N_COL; ++jcol){
      col = colIndex[jcol] - 1; // remove 1 since R indices (input) are 1-based
      for (int j = p[col]; j < p[col + 1]; ++j){
        double diff = x[j] - K;
        exv[i[j]] += diff;
        ex2v[i[j]] += diff*diff;
        count_zero_by_row[i[j]]--;
      }
    }

    // also remove 1 since R indices (input) are 1-based
    exv = exv[H0_genes-1];
    ex2v = ex2v[H0_genes-1];
    count_zero_by_row = count_zero_by_row[H0_genes-1];

    for(int j = 0; j < N_GENES; ++j){
      exv[j] += -1.0*K*count_zero_by_row[j];
      ex2v[j] += K*K*count_zero_by_row[j];
    }
    double ex = ( std::accumulate(exv.begin(), exv.end(), 0.0) ) ;
    double ex2 = ( std::accumulate(ex2v.begin(), ex2v.end(), 0.0) ) ;

    H0_sd = sqrt( (ex2 - (ex*ex) / n) / (n - 1.0) );
    H0_sd = H0_sd / sqrt(N_COL*1.0);

    // create a List containing the main results
    Rcpp::List results = Rcpp::List::create(
      Rcpp::Named("lig_avg_expr") = ligands_avg,
      Rcpp::Named("rec_avg_expr") = receptors_avg,
      Rcpp::Named("H0_mean") = H0_mean,
      Rcpp::Named("H0_sd") = H0_sd
    );

    return results;
  } // method BigBootstrap

  Rcpp::List BigBootstrapParallel(Rcpp::NumericVector colIndex,
                                  Rcpp::NumericVector H0_genes,
                                  Rcpp::NumericVector ligands,
                                  Rcpp::NumericVector receptors){

    const int N_COL = colIndex.length();
    const int N_GENES = H0_genes.length();
    const int N_LIG = ligands.length();
    const int N_REC = receptors.length();

    // compute the mean of all genes
    Rcpp::NumericVector geneMeans = rowMeansColIndex(colIndex);

    // compute ligands and receptors means (subset of all genes mean)
    Rcpp::NumericVector ligands_avg = geneMeans[ligands-1];
    Rcpp::NumericVector receptors_avg = geneMeans[receptors-1];

    // compute H0 genes means (subset of all genes mean) to be used in matrix mean computation
    Rcpp::NumericVector H0geneMeans = geneMeans[H0_genes-1];
    // compute matrix mean as the mean of genes means
    double H0_mean = ( std::accumulate(H0geneMeans.begin(), H0geneMeans.end(), 0.0) ) / N_GENES;


    double H0_sd = 0.0;
    double K = H0_mean; // shift parameter
    double n = N_GENES * N_COL; // total number of elements
    Rcpp::NumericVector exv(Dim[0]);
    Rcpp::NumericVector ex2v(Dim[0]);
    Rcpp::NumericVector count_zero_by_row(Dim[0], N_COL*1.0);

    int col;
    for (int jcol = 0; jcol < N_COL; ++jcol){
      col = colIndex[jcol] - 1; // remove 1 since R indices (input) are 1-based
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int j = p[col]; j < p[col + 1]; ++j){
        double diff = x[j] - K;
        exv[i[j]] += diff;
        ex2v[i[j]] += diff*diff;
        count_zero_by_row[i[j]]--;
      }
    }

    // also remove 1 since R indices (input) are 1-based
    exv = exv[H0_genes-1];
    ex2v = ex2v[H0_genes-1];
    count_zero_by_row = count_zero_by_row[H0_genes-1];

    for(int j = 0; j < N_GENES; ++j){
      exv[j] += -1.0*K*count_zero_by_row[j];
      ex2v[j] += K*K*count_zero_by_row[j];
    }
    double ex = ( std::accumulate(exv.begin(), exv.end(), 0.0) ) ;
    double ex2 = ( std::accumulate(ex2v.begin(), ex2v.end(), 0.0) ) ;

    H0_sd = sqrt( (ex2 - (ex*ex) / n) / (n - 1.0) );
    H0_sd = H0_sd / sqrt(N_COL*1.0);

    // create a List containing the main results
    Rcpp::List results = Rcpp::List::create(
      Rcpp::Named("lig_avg_expr") = ligands_avg,
      Rcpp::Named("rec_avg_expr") = receptors_avg,
      Rcpp::Named("H0_mean") = H0_mean,
      Rcpp::Named("H0_sd") = H0_sd
    );

    return results;
  } // method BigBootstrapParallel


  std::vector<std::vector<double>> extractCol(Rcpp::NumericVector colIndex) {
    std::vector<std::vector<double>> cols(Dim[0]);
    const int NCOL = colIndex.size();
    int col;
    for (int jcol = 0; jcol < NCOL; ++jcol){
      col = colIndex[jcol] - 1; // remove 1 since R indices (input) are 1-based
      for (int j = p[col]; j < p[col + 1]; ++j){
        cols[i[j]].push_back(x[j]);
      }
    }
    return cols;
  } // method extractCol
  
  
  
  
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ADDED METHODS
  
  
  
  void extractCol_and_col_idx(std::vector<std::vector<double>>& cols,
                              std::vector<std::vector<double>>& col_idx,
                              Rcpp::NumericVector colIndex) {
    //std::vector<std::vector<double>> cols(Dim[0]);
    const int NCOL = colIndex.size();
    int col;
    for (int jcol = 0; jcol < NCOL; ++jcol){
      col = colIndex[jcol] - 1; // remove 1 since R indices (input) are 1-based
      for (int j = p[col]; j < p[col + 1]; ++j){
        if(x[j] != 0){
          cols[i[j]].push_back(x[j]);
          col_idx[i[j]].push_back(col);
        }
      }
    }
    //return cols;
  } // method extractCol_and_col_idx
  
  
  
  
  // sign function. Returns the sign of the number given as input. Returns 0 if input = 0
  template <typename T> 
  int sgn(T val) {
    return (T(0) < val) - (val < T(0));
  }

  template <typename T>
  std::vector<size_t> sort_indexes(std::vector<T> const &v, bool abs = 0){
    
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    
    if(!abs){
      // sort indexes based on comparing values in v
      std::sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    }else{
      // sort indexes based on comparing values in abs(v). (abs = absolute value)
      std::sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return std::abs(v[i1]) < std::abs(v[i2]);});
    }
    
    
    return idx;
  }// method sort_indexes
  
  
  
  template <typename T>
  std::vector<double> calculate_ranks(const std::vector<T>& v, double& correctionFactor) {
    
    std::vector<double> ranks(v.size());
    
    size_t i=0, j=0;
    while (i < v.size()){
      j = i + 1;
      while (j < v.size()){
        if(v[i] != v[j]){
          break;
        }
        j++;
      }
      
      if((j-i) > 1){
        //countOfElements.push_back(static_cast<double>(j-i));
        double temp = j-i;
        correctionFactor += temp*temp*temp -temp;
      }
      
      for(size_t k = i; k <= j-1; k++){   
        ranks[k] = 1 + static_cast<double>(i + j-1)/static_cast<double>(2);
      }
      i = j;
    }
    
    
    return ranks;
  }// method calculate_ranks
  

  // sparse_diff returns V1-V2 taking into consideration their column indices in the sparse matrix. 
  // It omits zeroes in the result
  Rcpp::NumericVector sparse_diff(Rcpp::NumericVector V1,
                                  Rcpp::NumericVector V1_idx,
                                  Rcpp::NumericVector V2,
                                  Rcpp::NumericVector V2_idx){
    Rcpp::NumericVector X;
    
    
    // std::cout<<"\nV1: ";
    // for(size_t i=0; i<V1.size(); i++)
    //   std::cout<<V1[i]<<" ";
    // 
    // std::cout<<"\nV1_idx: ";
    // for(size_t i=0; i<V1_idx.size(); i++)
    //   std::cout<<V1_idx[i]<<" ";
    // 
    // std::cout<<"\nV2: ";
    // for(size_t i=0; i<V2.size(); i++)
    //   std::cout<<V2[i]<<" ";
    // 
    // std::cout<<"\nV2_idx: ";
    // for(size_t i=0; i<V2_idx.size(); i++)
    //   std::cout<<V2_idx[i]<<" ";
    
    size_t i=0, j=0;
    
    while(i<V1_idx.size() && j<V2_idx.size()){
      if(V1_idx[i] == V2_idx[j]){
        if((V1[i] - V2[j]) != 0)
          X.push_back(V1[i] - V2[j]);
        i++;
        j++;
      }else if(V1_idx[i]<V2_idx[j]){
        X.push_back(V1[i]);
        i++;
      }else{
        X.push_back(-V2[j]);
        j++;
      }
    }
    
    for(;i<V1_idx.size();i++)
        X.push_back(V1[i]);
      for(;j<V2_idx.size();j++)
          X.push_back(-V2[j]);
        
    
    
    
    
    return X;
  }// method sparse_diff
  
  
  
  void wilcox_no_zeroes (Rcpp::NumericVector V1,
                         Rcpp::NumericVector V2,
                         size_t n_zeroes_col_id,
                         size_t n_zeroes_ref_col_id,
                         double& bothtails,
                         double& righttail,
                         double& lefttail,
                         bool verbose = 0){
    
    std::vector <double> X;
    
    X.insert(X.end(), V1.begin(), V1.end());
    X.insert(X.end(), V2.begin(), V2.end());
    
    //total n. of zeroes
    double tot_0 = n_zeroes_col_id + n_zeroes_ref_col_id;
    
    int n = X.size() + tot_0;
    //int n = col_id.size() + ref_col_id.size();
    if(verbose)
      std::cout<<"\nn :"<<n;
    
    std::vector <double> ranks (X.size());
    
    //sort the vector and save the indices vector
    std::vector<size_t> idx = sort_indexes(X);
    
    sort(X.begin(), X.end());
    
    double correctionFactor = (tot_0*tot_0*tot_0)-tot_0;
    
    // fast version of getranks  
    ranks = calculate_ranks(X, correctionFactor);
    
    if(verbose){
      std::cout<<"\nranks: ";
      for(size_t i=0; i<ranks.size(); i++)
        std::cout<<ranks[i]<<" ";
    }
    
    
    
    
    //ranked sums
    long double T0 = 0, T1 = 0;
    
    for(size_t i=0; i<X.size(); i++){
      //idx[i]<V1.size()? T0 += ranks[i] : T1 += ranks[i];
      if(idx[i]<V1.size())
        T0 += ranks[i];
    }
    
    // aggiungo a T0 la somma dei rank dei primi n_zeroes_col_id numeri (cioè i ranks corrispondenti agli 0 di col_id) + 
    // + conto il fatto che i ranks trovati in precedenza devono essere traslati di n_zeroes_col_id
    
    
    //rank of each zero:
    double rank_0 = 0;
    if(tot_0 != 0)
      rank_0 = (  (tot_0*(tot_0+1))/2  )/tot_0;
    // std::cout<<"\nrank0: "<<rank_0;
    
    T0+=rank_0*n_zeroes_col_id + (tot_0)*(V1.size());
    
    T1 = (n * (n+1) / 2) - T0;
    
    
    double n0 = V1.size() +n_zeroes_col_id , n1 = V2.size() +n_zeroes_ref_col_id;
    
    if(verbose){
      std::cout<<"\nn0: "<<n0<<"\n";
      std::cout<<"n1: "<<n1<<"\n";
      std::cout<<"T0: "<<T0<<"\n";
      std::cout<<"T1: "<<T1<<"\n\n";
    }
    
    //Calculating U values:
    double long U0, U1;
    
    U0 = (n0*(n0+1)*(-0.5))+T0;
    U1 = (n1*(n1+1)*(-0.5))+T1;
    
    //double U; // U= min(U0, U1)
    //U0<U1 ? U = U0 : U = U1;
    if(verbose){
      std::cout<<"U0: "<<U0<<"\n\n";
      std::cout<<"U1: "<<U1<<"\n\n";
    }
    
    double muU = n0*n1/2;; //Expected value of U
    if(verbose)
      std::cout<<"muU: "<<muU<<"\n";
    
    //sigmaU = (sqrt(n0)*sqrt(n1)*sqrt((n0+n1+1)))/sqrt(12);  //NOT ADJUSTED FOR TIE CORRECTION!!!
    
    
    //standard error of U
    long double sigmaU = sqrt(((n0*n1)/12)*((n0+n1+1)-(  (correctionFactor)/((n0+n1)*(n0+n1-1))   )   ));
    if(verbose){
      std::cout<<"correctionFactor:"<<correctionFactor;
      std::cout<<"sigmaU: "<<sigmaU<<"\n";
    }
    
    //Z value
    //double z = (U-muU)/sigmaU;
    double z=0;
    z = U0<U1 ? U0-muU : U1-muU;
    if(z!=0)
      z = z<0 ? z+0.5 : z-0.5; // Continuity correction
    //std::cout<<"U-muU: "<<z<<"\n";
    if(sigmaU != 0)
      z = z/sigmaU;
    else
      z=0;
    if(verbose)
      std::cout<<"Z value: "<<z<<"\n";
    
    
    // returning the p values in the passed variables
    lefttail = (1+erf(z/sqrt(2)))/2;
    
    bothtails = 2*lefttail;  
    if(verbose)
      std::cout<<"bothtails: "<<bothtails<<"\n";
    
    righttail = 1-lefttail;
    
    // Rcpp::Rcout<<"P-value bothTails: "<<bothtails<<"\n";
    // Rcpp::Rcout<<"P-value rightTail: "<<righttail<<"\n";
    // Rcpp::Rcout<<"P-value leftTail : "<<lefttail<<"\n";
    // 
    
    
  }// method wilcox_no_zeroes
  
  
  // X contains the difference (omitting zeroes) between two vectors
  void wilcox_no_zeroes_paired(Rcpp::NumericVector X,
                               double& bothtails,
                               double& righttail,
                               double& lefttail,
                               bool verbose = 0){
    
    if(verbose){
      std::cout<<"\nX: ";
      for(size_t i=0; i<X.size(); i++)
        std::cout<<X[i]<<" ";
    }
    
    std::sort(X.begin(), X.end(), [](double i, double j) { return std::abs(i) < std::abs(j); });
    if(verbose){
      std::cout<<"\nsortedX: ";
      for(size_t i=0; i<X.size(); i++)
        std::cout<<X[i]<<" ";
    }
    
    std::vector<double> abs_X (X.size());
    for(size_t i=0; i<X.size(); i++)
      abs_X[i] = fabs(X[i]);
    
    if(verbose){
      std::cout<<"\nabs_X: ";
      for(size_t i=0; i<abs_X.size(); i++)
        std::cout<<abs_X[i]<<" ";
    }
    
    std::vector <double> ranks (X.size());
    double correctionFactor = 0;
    // fast version of getranks  
    ranks = calculate_ranks(abs_X, correctionFactor);
    
    if(verbose){
      std::cout<<"\nranks: ";
      for(size_t i=0; i<ranks.size(); i++)
        std::cout<<ranks[i]<<" ";
    }
    double T = 0;
    double Tplus = 0;
    double Tminus = 0;
    
    if(verbose)
      std::cout<<"\n";
    
    //calculate the T statistic
    for(size_t i=0; i<ranks.size(); i++){
      T += sgn(X[i])*ranks[i];
      sgn(X[i]) >0? Tplus += sgn(X[i])*ranks[i] : Tminus -= sgn(X[i])*ranks[i];
    //   if(verbose)
    //     std::cout<<sgn(X[i])*ranks[i]<<" ";
    }
    
    
    if(verbose){
      std::cout<<"\nT:  "<<T;
      std::cout<<"\nT+: "<<Tplus;
      std::cout<<"\nT-: "<<Tminus;
      std::cout<<"\n";
    }
    
    // Calculate sigma, mu and the p value
    
    //number of pairs with a non zero difference
    double n = X.size();
    
    //long double sigma = sqrt(  ((n*(n+1)*((2*n)+1)) - (nz*(nz+1)*(2*nz+1)) - (correctionFactor/2))  /24);
    long double sigma = sqrt(  ((n*(n+1)*((2*n)+1)) - (correctionFactor/2) )  /24);
    if(verbose)
      std::cout<<"\nsigma: "<<sigma;
    // std::cout<<"\nzero correction: "<<(nz*(nz+1)*(2*nz+1));
    // std::cout<<"\nties correction: "<<(correctionFactor/2);  
    
    //double mu = ((n*(n+1))/4) - ((nz*(nz+1))/4);
    double mu = ((n*(n+1))/4);
    if(verbose){
      std::cout<<"\nmu: "<<mu;
      std::cout<<"\nn: "<<n;
    }
    double z = Tplus<Tminus ? Tplus-mu : Tminus-mu;
    if(z!=0)
      z = z<0 ? z+0.5 : z-0.5; // Continuity correction
    if(sigma != 0)
      z = z/sigma;
    else
      z = 0;
    if(verbose)
      std::cout<<"\nz: "<<z<<"\n";
    //std::cout<<"Z value: "<<z<<"\n";
    
    
    // returning the p values in the passed variables
    lefttail = (1+erf(z/sqrt(2)))/2;
    
    righttail = 1-lefttail;
    
    bothtails = 2*((lefttail < righttail) ? lefttail : righttail); 
    if(verbose)
      std::cout<<"\nbothtails: "<<bothtails<<"\n";
    ////////////////////////////////////////////////////////////////
    
    
   
    if(verbose)
      std::cout<<"\n*** END of sparseWilcox_paired_single_row ***\n";
    
  }// method wilcox_no_zeroes_paired
  

  Rcpp::DataFrame sparseWilcox (Rcpp::NumericVector rowIndex,
                                           Rcpp::NumericVector col_id, 
                                           Rcpp::NumericVector ref_col_id,
                                           bool paired,
                                           std::string alternative = "two.sided",
                                           bool verbose = 0){
    
    if(verbose)
      std::cout<<"\n*** sparseWilcox ***\n";
    
    Rcpp::NumericVector pvalues(nrow(), 1.0);
    Rcpp::NumericVector sign(nrow());
    sign.fill(NA_REAL); // NA_REAL is a constant representing NaN in Rcpp
    
    
    double bothtails = -3;
    double righttail = -3;
    double lefttail = -3;
    Rcpp::NumericVector pValues(rowIndex.size());
    
    

    

      // std::vector<std::vector<double>>  tmp1s = extractCol(col_id);
      // std::vector<std::vector<double>> tmp2s = extractCol(ref_col_id);

    
    std::vector<std::vector<double>>  tmp1s(Dim[0]);
    std::vector<std::vector<double>>  tmp1s_idx(Dim[0]);
    extractCol_and_col_idx(tmp1s, tmp1s_idx, col_id);
    
    std::vector<std::vector<double>>  tmp2s(Dim[0]);
    std::vector<std::vector<double>>  tmp2s_idx(Dim[0]);
    extractCol_and_col_idx(tmp2s, tmp2s_idx, ref_col_id);

    
    for (int rr = 0; rr < rowIndex.length(); rr++) {
      
      int ii = rowIndex[rr] -1 ;
      
      Rcpp::NumericVector tmp1 = Rcpp::wrap(tmp1s[ii]);
      Rcpp::NumericVector tmp1_idx = Rcpp::wrap(tmp1s_idx[ii]);
      
      // std::cout<<"\ntmp1: \n";
      // for(size_t i=0; i<tmp1.size(); i++)
      //   std::cout<<tmp1[i]<<" ";
      // 
      // std::cout<<"\ntmp1_idx: \n";
      // for(size_t i=0; i<tmp1_idx.size(); i++)
      //   std::cout<<tmp1_idx[i]<<" ";
      
      //adding zeroes
      // int z = col_id.size() - tmp1.size();
      // for (int k=0; k < z;  ++k){
      //   tmp1.push_back(0.0);
      // }
      
      double sum1 = std::accumulate(tmp1.begin(), tmp1.end(), 0.0);
      double mean1 = sum1 / col_id.size();
      
      Rcpp::NumericVector tmp2 = Rcpp::wrap(tmp2s[ii]);
      Rcpp::NumericVector tmp2_idx = Rcpp::wrap(tmp2s_idx[ii]);
      
      size_t temp = ref_col_id[0]-1;
      for(size_t i=0; i<tmp2_idx.size(); i++)
        tmp2_idx[i] -= temp;
      
      
      //adding zeroes
      // z = ref_col_id.size() - tmp2.size();
      // for (int k=0; k < z;  ++k){
      //   tmp2.push_back(0.0);
      // }
      
      // NumericVector tmp2 = row[ref_col_id -1];
      double sum2 = std::accumulate(tmp2.begin(), tmp2.end(), 0.0);
      double mean2 = sum2 / ref_col_id.size();
      
      if(mean2 <= mean1)
      {
        sign[ii] = 1;
      }else{
        sign[ii] = 0;
      }
      
      size_t n_zeroes_col_id = col_id.size() - tmp1.size();
      size_t n_zeroes_ref_col_id = ref_col_id.size() - tmp2.size();
      
      if(!paired){
        wilcox_no_zeroes(tmp1,
                         tmp2,
                         n_zeroes_col_id,
                         n_zeroes_ref_col_id,
                         bothtails,
                         righttail,
                         lefttail,
                         verbose);
        }
      else{
        if(col_id.length() != ref_col_id.length())
          throw(Rcpp::exception("'col_id' and 'ref_col_id' must have the same length"));
        else{
          Rcpp::NumericVector X = sparse_diff(tmp1, tmp1_idx, tmp2, tmp2_idx);
          
          wilcox_no_zeroes_paired(X,
                                  bothtails,
                                  righttail,
                                  lefttail,
                                  verbose);
        }
      }
      
      if(alternative == "two.sided"){
        pvalues[ii] = bothtails;
        //std::cout<<bothtails;
        }
      if(alternative == "less")
        pvalues[ii] = righttail;
      if(alternative == "greater")
        pvalues[ii] = lefttail;
      
    }
    
    
      
      // Creating DataFrame to return
      Rcpp::DataFrame results = Rcpp::DataFrame::create( Rcpp::Named("p.value") = pvalues,
                                                   Rcpp::Named("sign") = sign);
    
    
    if(verbose)
      std::cout<<"\n*** END of sparseWilcox ***\n";
    
    return results;
    
  }// method sparseWilcox
  
  //END OF ADDED METHODS
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Rcpp::NumericVector extractRow(int row){
    // row = row- 1;

    // Rcpp::NumericVector r(Dim[1], 0.0);
    //
    // auto it = std::find(i.begin(), i.end(), row);
    //
    // for (int col = 0; col < Dim[1]; ++col) {
    //   if(p[col] < std::distance(i.begin(), it)){
    //     continue;
    //   }
    //
    //   for (int j = p[col]; j < p[col + 1]; ++j) {
    //     if (i[j] == row) r[col] = x[j];
    //     else if (i[j] > row) break;
    //   }
    //
    //   it = std::find(it + 1, i.end(), row);
    //   if(it == i.end()){
    //     break;
    //   }
    // }

    Rcpp::NumericVector r(Dim[1], 0.0);



    for (int col = 0; col < Dim[1]; ++col) {
      for (int j = p[col]; j < p[col + 1]; ++j) {
        if (i[j] == row) r[col] = x[j];
        else if (i[j] > row) break;
      }

    }


    // const int NCOL = colIndex.size();
    // int col;
    // for (int jcol = 0; jcol < NCOL; ++jcol){
    //    col = colIndex[jcol]; // remove 1 since R indices (input) are 1-based
    //   for (int j = p[jcol]; j < p[jcol + 1]; ++j){
    //     if (i[j] == row) r[jcol] = x[j];
    //     else if (i[j] > row) break;
    //   }
    // }


    return r;
  }



};
//}

