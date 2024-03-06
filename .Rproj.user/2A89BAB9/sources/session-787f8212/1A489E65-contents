//AUTHOR: Lorenzo Benatti
//e-mail: lorenzo.benatti98@gmail.com



#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]]



// Other functions used
//###################################################################################################################

// Simple function to display a 1-d vector
template <typename T>
void display(std::vector<T> const &vector){

  for (int i = 0; i < vector.size(); i++)
    std::cout <<std::fixed << std::setprecision(2) << std::setw(5) << std::setfill(' ') << vector[i] << ' ';

  std::cout<<"\n";
}

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
}




// takes as input a sorted vector
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
}







void mannwhitneyutest( 
    Rcpp::NumericVector const &V1,
    Rcpp::NumericVector const &V2,
    double& bothtails,
    double& righttail,
    double& lefttail
)
{
  int V1_size = V1.size();
  //copy V1 and V2 into X
  std::vector <double> X(V1.size());
  std::copy(V1.begin(), V1.end(), X.begin());
  X.insert( X.end(),V2.begin(),V2.end());
  
  
  
  int n = X.size();
  
  
  std::vector <double> ranks (X.size());
  //std::vector <double> countOfElements;   // This vector will store the number of ties for the K-th rank. (size(countOfElements) <  size(ranks))
  // we use countOfElements to account for tie correction when calculating sigmaU
  
  //sort the vector and save the indices vector
  std::vector<size_t> idx = sort_indexes(X);
  
  //for(size_t i=0;i< idx.size(); i++)
  //  std::cout<<idx[i]<<" ";
  
  // std::vector<double> sorted_copy(X.size());
  // double temp=0;
  // for(size_t i=0; i<X.size(); i++){
  //  sorted_copy[i] = X[idx[i]];
  // }
  sort(X.begin(), X.end());
  
  double correctionFactor = 0;
  // fast version of getranks  
  ranks = calculate_ranks(X, correctionFactor);
  
  
  // std::cout<<"ranks : ";
  // display(ranks);
  
  // std::cout<<"countOfElements : ";
  // display(countOfElements);
  
  
  
  //ranked sums
  long double T0 = 0, T1 = 0;
  
  for(size_t i=0; i<X.size(); i++){
    //idx[i]<V1.size()? T0 += ranks[i] : T1 += ranks[i];
    if(idx[i]<V1.size())
      T0 += ranks[i];
  }
  
  T1 = (X.size() * (X.size()+1) / 2) - T0;
  
  
  double n0 = V1.size(), n1 = V2.size();
  
  //now i can restore X to the original size
  //X.erase(X.end()-n1, X.end());
  
  // std::cout<<"n0: "<<n0<<"\n";
  // std::cout<<"n1: "<<n1<<"\n";
  //std::cout<<"T0: "<<T0<<"\n";
  //std::cout<<"T1: "<<T1<<"\n\n";
  
  //Calculating U values:
  double long U0, U1;
  
  //U0 = n0*n1 + (n0*(n0+1)/(2))-T0;    //Per wikipedia U0 = T0-((n0*(n0+1)/2))
  //U1 = n0*n1 + (n1*(n1+1)/(2))-T1;    //Per wikipedia U1 = T1-((n1*(n1+1)/2))
  //n0*n1 = U0+U1 -> sempre detto da wikipedia
  
  U0 = (n0*(n0+1)*(-0.5))+T0;
  U1 = (n1*(n1+1)*(-0.5))+T1;
  
  //double U; // U= min(U0, U1)
  //U0<U1 ? U = U0 : U = U1;
  //std::cout<<"U0: "<<U0<<"\n\n";
  //std::cout<<"U1: "<<U1<<"\n\n";
  
  double muU = n0*n1/2;; //Expected value of U
  //std::cout<<"muU: "<<muU<<"\n";
  
  
  
  //sigmaU = (sqrt(n0)*sqrt(n1)*sqrt((n0+n1+1)))/sqrt(12);  //NOT ADJUSTED FOR TIE CORRECTION!!!
  
  
  //double correctionFactor = 0;
  //for (int i=0; i<countOfElements.size(); i++){
  //  correctionFactor += pow(countOfElements[i],3) - countOfElements[i];
  //}
  
  
  //if(!countOfElements.empty()){
  //  for(std::vector<double>::iterator it = countOfElements.begin(); it != countOfElements.end(); ++it){
  //    correctionFactor += ( ((*it) * (*it) * (*it)) - *it );
  //  }
  //}
  
  
  // cout<<"n. of ties:"<<size(countOfElements)<<"\n";
  // display(countOfElements)
  
  //standard error of U
  long double sigmaU = sqrt(((n0*n1)/12)*((n0+n1+1)-(  (correctionFactor)/((n0+n1)*(n0+n1-1))   )   ));
  //std::cout<<"sigmaU: "<<sigmaU<<"\n";
  
  
  //Z value
  //double z = (U-muU)/sigmaU;
  double z=0;
  z = U0<U1 ? U0-muU : U1-muU;
  z = z<0 ? z+0.5 : z-0.5; // Continuity correction
  //std::cout<<"U-muU: "<<z<<"\n";
  z = z/sigmaU;
  //std::cout<<"Z value: "<<z<<"\n";
  
  
  // returning the p values in the passed variables
  lefttail = (1+erf(z/sqrt(2)))/2;
  
  bothtails = 2*lefttail;  
  
  righttail = 1-lefttail;
  
}

void wilcoxonSignedRankTest( 
    Rcpp::NumericVector const &V1,
    Rcpp::NumericVector const &V2,
    double& bothtails,
    double& righttail,
    double& lefttail,
    bool verbose = 0
){
  std::vector <double> X;

  //numero di zeri
  size_t nz = 0;

  for (size_t i=0; i<V1.size(); i++){
    if((V1[i]- V2[i]) != 0)
      X.push_back(V1[i]- V2[i]);
    else
      nz++;
    }
  
    if(verbose){
      for(size_t i=0; i<X.size(); i++)
      std::cout<<X[i]<<" ";
    
      std::cout<<"\nnz: "<<nz<<"\n";
      }
  

  sort(X.begin(), X.end(), [](double i, double j) { return std::abs(i) < std::abs(j); });
  if(verbose){
    for(size_t i=0; i<X.size(); i++)
      std::cout<<X[i]<<" ";
  }
  

  
  std::vector <double> ranks (X.size());
  double correctionFactor = 0;
  // fast version of getranks  
  ranks = calculate_ranks(X, correctionFactor);
  
  double T = 0;
  double Tplus = 0;
  double Tminus = 0;
  
  if(verbose)
    std::cout<<"\n";
  
  //calculate the T statistic
  for(size_t i=0; i<ranks.size(); i++){
    T += sgn(X[i])*ranks[i];
    sgn(X[i]) >0? Tplus += sgn(X[i])*ranks[i] : Tminus -= sgn(X[i])*ranks[i];
    if(verbose)
      std::cout<<sgn(X[i])*ranks[i]<<" ";
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
    z=0;
  if(verbose)
  std::cout<<"\nz: "<<z<<"\n";
  //std::cout<<"Z value: "<<z<<"\n";
  
  
  // returning the p values in the passed variables
  lefttail = (1+erf(z/sqrt(2)))/2;
  
  righttail = 1-lefttail;
  
  bothtails = 2*(lefttail < righttail ? lefttail : righttail); 
  
}




//#############################################################################################################
// This is the actual function visible by R

// [[Rcpp::export]]
DataFrame fasterWilcox_v3(NumericVector V1, NumericVector V2, bool verbose = 0, const std::string& alternative = "two.sided", bool paired = 0) {
  
  
  double bothTails = -3;
  double rightTail = -3;
  double leftTail  = -3;
  
  double results = -3;
  
  if(!paired)
    mannwhitneyutest(V1, V2, bothTails, rightTail, leftTail);
  else{
    if(V1.length() != V2.length())
      throw(Rcpp::exception("'x' and 'y' must have the same length"));
    else
      wilcoxonSignedRankTest(V1, V2, bothTails, rightTail, leftTail, verbose);
    }
  
  //results[0] = bothTails;
  //results[1] = leftTail;
  //results[2] = rightTail;
  
  if(verbose){
    Rcout<<"\n";
    Rcout<<"Mann Whitney U test results: \n";
    Rcout<<"P-value bothTails: "<<bothTails<<"\n";
    Rcout<<"P-value rightTail: "<<rightTail<<"\n";
    Rcout<<"P-value leftTail : "<<leftTail<<"\n";
  }
  
  
  // Creating DataFrame df
  //DataFrame res = DataFrame::create( Named("p.value") = bothTails,         
  //                                  Named("p.leftTail") = rightTail,
  //                                  Named("p.rightTail") = leftTail); 
  
  if(alternative == "two.sided")
    results = bothTails;
  if(alternative == "greater")
    results = leftTail;
  if(alternative == "less")
    results =  rightTail;
  
  // Creating DataFrame df
  DataFrame res = DataFrame::create( Named("p.value") = results);  
  
  return res;
  
}










/*** R
#rcpp_MWU_Lor(c(45, 33, 35, 39, 42), c(34, 36, 41, 43, 44, 37))
*/