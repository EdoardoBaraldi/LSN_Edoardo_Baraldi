#include "mylib.h"

using namespace std;

int main(){

 vector<double> sum_prog(N,0.);
 vector<double> sum2_prog(N,0.);
 vector<double> err_prog(N,0.);
 double sum=0.;	
 
 	
 Random Rand;			//random generator
 pray_the_sun(Rand); 		//initialize random generator
 //ex 1.1	
 for(int i=0; i<N; i++){				//cycle on N blocks	
  for(int j=0; j<L; j++){				//cycle on the L elements of each block
   sum +=Rand.Rannyu();}				//sum of L elements
  if(i==0){
   sum_prog[i]=sum/double(L);				//mean of block 0
   sum2_prog[i]=pow(sum/double(L),2);			//square mean of block 0
  }
  else{
   sum_prog[i]=sum_prog[i-1] + sum/double(L);  	        // sum of r_j with {j=0,...,i}
   sum2_prog[i]=sum2_prog[i-1] + pow(sum/double(L),2);  // sum of (r_j)^2 with {j=0,...,i}
  }
  sum=0;         					//reset sum to 0 for new block	
 }
  for(int i=0;i<N; i++){
   sum_prog[i] /= (i+1);				//progressive mean of blocks
   sum2_prog[i] /= (i+1);				//progressive square mean of blocks
   err_prog[i] = error (sum_prog, sum2_prog, i);	//progressive error of blocks
  }	
  
  Printfile_2( sum_prog, err_prog, "Data1_1");		//print in a file N progressive mean with their own error
  
  //ex 1.2
  
  for(int i=0; i<N; i++){				//cycle on N blocks	
   for(int j=0; j<L; j++){				//cycle on the L elements of each block
    sum +=pow(Rand.Rannyu()-0.5,2);}			//sum of (random numbers - 0.5)^2
   if(i==0){
    sum_prog[i]=sum/double(L);				//mean of block 0
    sum2_prog[i]=pow(sum/double(L),2);			//square mean of block 0
   }
   else{
    sum_prog[i]=sum_prog[i-1] + sum/double(L);  	// sum of r_j with {j=0,...,i}
    sum2_prog[i]=sum2_prog[i-1] + pow(sum/double(L),2);  // sum of (r_j)^2 with {j=0,...,i}
   }
  sum=0;         					//reset sum to 0 for new block	
  }
  for(int i=0;i<N; i++){
   sum_prog[i] /= (i+1);				//progressive mean of blocks
   sum2_prog[i] /= (i+1);				//progressice square mean of blocks
   err_prog[i] = error (sum_prog, sum2_prog, i);	//progressive error of blocks
  }	
  
  Printfile_2( sum_prog, err_prog, "Data1_2");		//print in a file N progressive mean with their own error
  
  //ex 1.3
  
  int Mchi=100;
  int nchi=pow(10,4);
  double events=nchi/Mchi;
  int times=100;
  vector<double> chi_squared(times,0.);
  
  for(int i=0; i<times; i++){
   vector<double> counts(Mchi,0.);
   for(int j=0; j<nchi; j++){
    counts[floor(Rand.Rannyu()*100)] ++;}
   for(int k=0;k<Mchi; k++){
    chi_squared[i] += pow(counts[k] - events, 2)/events;}
  }
  
  Printfile(chi_squared, "Data1_3");

return 0;	
}
