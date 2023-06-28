#include "mylib.h"

using namespace std;

const double T = 1.;	// delivery time
const double K = 100.;	// strike price
const double r = 0.1;	// risk-free interest rate
const double volat = 0.25;	//volatility
const double S0 = 100.;	// asset price at t=0
const double t0 = 0;
const int step = 100;	//numebr of steps for discretized

double max(const double a, const double b) {
	return (a > b) ? a : b;
}

int main(){

 vector<double> sum_C(N,0.);
 vector<double> sum2_C(N,0.);
 vector<double> err_C(N,0.);
 
 vector<double> sum_P(N,0.);
 vector<double> sum2_P(N,0.);
 vector<double> err_P(N,0.);
 
 double sum_Call=0.;
 double sum_Put=0.;	
 
 	
 Random Rand;			//random generator
 pray_the_sun(Rand); 		//initialize random generator
 
 //ex 3.1	
 for(int i=0; i<N; i++){				//cycle on N blocks	
  for(int j=0; j<L; j++){				//cycle on the L elements of each block
   double S = S0*exp((r-pow(volat,2)/2)*T+volat*Rand.Gauss(0.,T));
   sum_Call += exp(-r*T)*max(0,S-K);			
   sum_Put += exp(-r*T)*max(0,K-S);
  }
  if(i==0){
   sum_C[i]=sum_Call/double(L);				//mean of block 0
   sum2_C[i]=pow(sum_Call/double(L),2);			//square mean of block 0
   sum_P[i]=sum_Put/double(L);				//mean of block 0
   sum2_P[i]=pow(sum_Put/double(L),2);
  }
  else{
   sum_C[i]=sum_C[i-1] + sum_Call/double(L);  	        // sum of call_j with {j=0,...,i}
   sum2_C[i]=sum2_C[i-1] + pow(sum_Call/double(L),2);  	// sum of call_j^2 with {j=0,...,i}
   sum_P[i]=sum_P[i-1] + sum_Put/double(L);  	        // sum of put_j with {j=0,...,i}
   sum2_P[i]=sum2_P[i-1] + pow(sum_Put/double(L),2);  	// sum of put_j^2 with {j=0,...,i}
  }
  sum_Call=0;
  sum_Put=0;        					//reset sum to 0 for new block	
 }
  for(int i=0;i<N; i++){
   sum_C[i] /= (i+1);				//progressive mean of blocks
   sum2_C[i] /= (i+1);				//progressive square mean of blocks
   err_C[i] = error (sum_C, sum2_C, i);		//progressive error of blocks
   sum_P[i] /= (i+1);				//progressive mean of blocks
   sum2_P[i] /= (i+1);				//progressive square mean of blocks
   err_P[i] = error (sum_P, sum2_P, i);		//progressive error of blocks 
  }	
  
  Printfile_2( sum_C, err_C, "Data3_1_C"); //print in a file N progressive mean with their own error
  Printfile_2( sum_P, err_P, "Data3_1_P");
 
 fill(sum_C.begin(), sum_C.end(), 0.);
 fill(sum2_C.begin(), sum2_C.end(), 0.);
 fill(sum_P.begin(), sum_P.end(), 0.);
 fill(sum2_P.begin(), sum2_P.end(), 0.);
 
 //ex 3.2
 
 double S_disc = S0;
 
 for(int i=0; i<N; i++){					
  for(int j=0; j<L; j++){
   for(int t=0; t<step; t++){			
    S_disc *= exp((r-pow(volat,2)/2)*(T/double(step))+volat*Rand.Gauss(0.,1.)*sqrt(T/double(step)));}
   cout << S_disc << endl;
   sum_Call += exp(-r*T)*max(0,S_disc-K);			
   sum_Put += exp(-r*T)*max(0,K-S_disc);
   S_disc = S0;
  }
  if(i==0){
   sum_C[i]=sum_Call/double(L);				//mean of block 0
   sum2_C[i]=pow(sum_Call/double(L),2);			//square mean of block 0
   sum_P[i]=sum_Put/double(L);				//mean of block 0
   sum2_P[i]=pow(sum_Put/double(L),2);
  }
  else{
   sum_C[i]=sum_C[i-1] + sum_Call/double(L);  	        // sum of call_j with {j=0,...,i}
   sum2_C[i]=sum2_C[i-1] + pow(sum_Call/double(L),2);  	// sum of call_j^2 with {j=0,...,i}
   sum_P[i]=sum_P[i-1] + sum_Put/double(L);  	        // sum of put_j with {j=0,...,i}
   sum2_P[i]=sum2_P[i-1] + pow(sum_Put/double(L),2);  	// sum of put_j^2 with {j=0,...,i}
  }
  sum_Call=0;
  sum_Put=0;        					//reset sum to 0 for new block	
 }
  for(int i=0;i<N; i++){
   sum_C[i] /= (i+1);				//progressive mean of blocks
   sum2_C[i] /= (i+1);				//progressive square mean of blocks
   err_C[i] = error (sum_C, sum2_C, i);		//progressive error of blocks
   sum_P[i] /= (i+1);				//progressive mean of blocks
   sum2_P[i] /= (i+1);				//progressive square mean of blocks
   err_P[i] = error (sum_P, sum2_P, i);		//progressive error of blocks 
  }	
  
  Printfile_2( sum_C, err_C, "Data3_1_C_disc"); //print in a file N progressive mean with their own error
  Printfile_2( sum_P, err_P, "Data3_1_P_disc");
 
 fill(sum_C.begin(), sum_C.end(), 0.);
 fill(sum2_C.begin(), sum2_C.end(), 0.);
 fill(sum_P.begin(), sum_P.end(), 0.);
 fill(sum2_P.begin(), sum2_P.end(), 0.);
 
return 0;

}
 
 
 
 
 
 
 
 
 
 
 
 
 
