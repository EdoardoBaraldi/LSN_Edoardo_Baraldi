#include "mylib.h"

using namespace std;

int main(){

 vector<double> sum_prog(N,0.);
 vector<double> sum2_prog(N,0.);
 vector<double> err_prog(N,0.);
 double sum=0.;	
 
 	
 Random Rand;			//random generator
 pray_the_sun(Rand); 		//initialize random generator
 
 //ex 2.1.1	
 for(int i=0; i<N; i++){				//cycle on N blocks	
  for(int j=0; j<L; j++){				//cycle on the L elements of each block
   sum +=(M_PI/2.*cos(M_PI*Rand.Rannyu()/2.));}		//sum of L elements of the function integrated
  if(i==0){
   sum_prog[i]=sum/double(L);				//mean of block 0
   sum2_prog[i]=pow(sum/double(L),2);			//square mean of block 0
  }
  else{
   sum_prog[i]=sum_prog[i-1] + sum/double(L);  	        // sum 
   sum2_prog[i]=sum2_prog[i-1] + pow(sum/double(L),2);  // squared sum
  }
  sum=0;         					//reset sum to 0 for new block	
 }
  for(int i=0;i<N; i++){
   sum_prog[i] /= (i+1);				//progressive mean of blocks
   sum2_prog[i] /= (i+1);				//progressive square mean of blocks
   err_prog[i] = error (sum_prog, sum2_prog, i);	//progressive error of blocks
  }	
  
  Printfile_2( sum_prog, err_prog, "Data1_1");		//print in a file N progressive mean with their own error
 
 fill(sum_prog.begin(), sum_prog.end(), 0.);
 fill(sum2_prog.begin(), sum2_prog.end(), 0.);
 //ex 2.1.2
 //i choose as a probability density a normalized straight line f(x)=2(1-x)
 
 for(int i=0; i<N; i++){				//cycle on N blocks	
  for(int j=0; j<L; j++){				//cycle on the L elements of each block
   double s=(1.-sqrt(1.-Rand.Rannyu()));	                        
   sum +=(M_PI/2.*cos(M_PI*s/2.)/(2.-2*s));}	//sum of L elements of the function
  if(i==0){
   sum_prog[i]=sum/double(L);				//mean of block 0
   sum2_prog[i]=pow(sum/double(L),2);			//square mean of block 0
  }
  else{
   sum_prog[i]=sum_prog[i-1] + sum/double(L);  	        // sum 
   sum2_prog[i]=sum2_prog[i-1] + pow(sum/double(L),2);  // squared sum
  }
  sum=0;         					//reset sum to 0 for new block	
 }
  for(int i=0;i<N; i++){
   sum_prog[i] /= (i+1);				//progressive mean of blocks
   sum2_prog[i] /= (i+1);				//progressive square mean of blocks
   err_prog[i] = error (sum_prog, sum2_prog, i);	//progressive error of blocks
  }
 
 Printfile_2( sum_prog, err_prog, "Data1_2");
 
 //ex 2.1.3
 //i choose as a probability density a normalized straight line f(x)=3x^2
 
 for(int i=0; i<N; i++){				//cycle on N blocks	
  for(int j=0; j<L; j++){				//cycle on the L elements of each block
   double s=cbrt(Rand.Rannyu());	                        
   sum +=(M_PI/2.*cos(M_PI*s/2.)/(3*s*s));}		//sum of L elements of the function
  if(i==0){
   sum_prog[i]=sum/double(L);				//mean of block 0
   sum2_prog[i]=pow(sum/double(L),2);			//square mean of block 0
  }
  else{
   sum_prog[i]=sum_prog[i-1] + sum/double(L);  	        // sum 
   sum2_prog[i]=sum2_prog[i-1] + pow(sum/double(L),2);  // squared sum
  }
  sum=0;         					//reset sum to 0 for new block	
 }
  for(int i=0;i<N; i++){
   sum_prog[i] /= (i+1);				//progressive mean of blocks
   sum2_prog[i] /= (i+1);				//progressive square mean of blocks
   err_prog[i] = error (sum_prog, sum2_prog, i);	//progressive error of blocks
  }
 
 Printfile_2( sum_prog, err_prog, "Data1_3");
 
return 0;

}
