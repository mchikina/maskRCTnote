#include <iostream>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <random>
#include <fstream>
#include <vector>
//#define NDEBUG
#include <cassert>
#include <Rcpp.h>
#include <RcppCommon.h>
#include <utility>		// USES pair

namespace Rcpp {
	template <typename T1, typename T2> SEXP wrap( const std::pair<T1,T2>& );
}

template<typename T1, typename T2> SEXP Rcpp::wrap( const
std::pair<T1,T2>& _v ) {
    return Rcpp::List::create(
        Rcpp::Named("first")  = Rcpp::wrap<T1>( _v.first ),
        Rcpp::Named("second") = Rcpp::wrap<T2>( _v.second )
	);
};

using namespace std;






void Tokenize( const char* szString, vector<std::string>&
               vecstrTokens, const char* szSeparators, bool fNoEmpties){
 const char* pc;
 string strCur;
 bool fPush;

 if( !( pc = szString ) )
  return;

 fPush = false;
 while( true ) {
   strCur.clear( );
   if( fNoEmpties )
    for( ; *pc && strchr( szSeparators, *pc ); ++pc );
   if( !*pc ) {
     if( !fNoEmpties && fPush )
      vecstrTokens.push_back( strCur );
     return; }
   for( ; *pc && !strchr( szSeparators, *pc ); ++pc )
    strCur += *pc;
   if( (fPush = !!*pc) )
    pc++;
   vecstrTokens.push_back( strCur );
  }
}


template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

// [[Rcpp::export]]
pair<vector<int>,double> permute(vector<double> A, vector<double> B, double min_effect, double max_effect, int bin_num, int T){
  double effect_range=max_effect-min_effect;

  vector<int> output(bin_num,0);

  int n=A.size();
  if (n!=B.size()){
    cerr << "Error: A and B have different sizes"<< endl;
    exit(-1);
  }

  mt19937_64 gen(314159265358979);
  bernoulli_distribution coin_flip(0.5);

  double initial_effect;
  {
    double total_A=0;
    double total_B=0;
    for (int i=0; i<n; i++){
      total_A+=A[i];
      total_B+=B[i];
    }
    initial_effect=(total_B-total_A)/total_B;
  }
  int num_total=0;
  int num_better=0;
  for (int iter=0; iter<T; iter++){
    double total_A=0;
    double total_B=0;
    for (int i=0; i<n; i++){
      if (coin_flip(gen)){
        total_A+=B[i];
        total_B+=A[i];
      }
      else{
        total_A+=A[i];
        total_B+=B[i];
      }
    }
    double effect=(total_B-total_A)/total_A;
    if (effect<min_effect)
      output[0]++;
    else if(effect>max_effect)
      output[bin_num-1]++;
    else
      output[1+(bin_num-2)*( (effect-min_effect)/effect_range)]++;
  //only use for p-value if the sign matches
 if(sgn(effect)!=sgn(initial_effect))
   continue;
 num_total++;
    if (abs(effect)>abs(initial_effect))
     num_better++;
    
  }
  return pair<vector<int>,double>(output,(double)num_better/(double)num_total);
}

int main () {
  string line;
  ifstream myfile ("input.txt");

  if (!myfile.good()){
    cerr << "ERROR with file" << endl;
    exit(-1);
  }
  getline(myfile, line);
  int n=atoi(line.c_str());

  vector<double> A(n);
  vector<double> B(n);

  int line_num=0;

  while ( getline (myfile,line) ){
    if (line_num==n){
     cout << "Ending after "<<n<< "lines";
     break;
    }
    vector<string> dataline;
    Tokenize(line.c_str(), dataline, "\t", true);
    if (dataline.size()!=2){
        cerr << "Error: expected row of length 2 in line "<< line_num<<", got length "<<dataline.size()<<endl;
        return -1;
    }
    A[line_num]=atof(dataline[0].c_str());
    B[line_num]=atof(dataline[1].c_str());
    line_num++;
  }

  int bins=100;
  int T=1000000;
  int barsize=100;
  pair<vector<int>,int> output=permute(A,B,-.5,.5,bins,T);
  assert(output.first.size()==bins);

  int output_max=0;
  for (int b=0; b<bins; b++)
   output_max=max(output_max,output.first[b]);


  for (int b=0; b<bins; b++){
    for (int i=0; i<barsize*output.first[b]/output_max; i++)
     cout << "*";
    cout << endl;
  }
  
  cout << "num_better is "<<output.second<<endl;
  cout << "p="<<(double)output.second/T<<endl;

  myfile.close();
  return 0;
}
