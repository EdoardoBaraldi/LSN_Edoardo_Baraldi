#include "mylib.h"

using namespace std;

const double Len=1.2;             	//lenght L of the stick
const double d=2.0;			//distance d between two straight lines, the number of throws is L

int main(){
  
 // i create vectors of value of pi, pi squared and error of pi with N dimensions
 vector<double> pi, pi2, err_pi;       
 pi.reserve(N);
 pi2.reserve(N);
 err_pi.reserve(N);
  
 double sum_hits = 0.0;
 double inv_LD = 2.0 * Len *L / d;
	
 
 Random Rand;			//random generator
 pray_the_sun(Rand); 		//initialize random generator

 for(int i=0; i<N; i++){	//cycle on N blocks	
  for(int j=0; j<L; j++){	//cycle on the L elements of each block
   if(count_hits(Rand, Len, d)){
    sum_hits+=1.0;}
  }
  if(i==0){
   pi.push_back(inv_LD/sum_hits);			//mean of block 0
   pi2.push_back(pi[i] * pi[i]);			//square mean of block 0
  }
  else{
   pi.push_back(pi[i - 1] + inv_LD/sum_hits); 	        // sum of p_j with {j=0,...,i}
   pi2.push_back(pi2[i - 1] + (inv_LD/sum_hits) * (inv_LD/sum_hits));   	// sum of (p_j)^2 with {j=0,...,i}
  }
  sum_hits=0.; 				
 }
 
 for(int i=0;i<N; i++){
   pi[i] /= (i+1);				//progressive mean of blocks
   pi2[i] /=  (i+1);			//progressive square mean of blocks
   err_pi.push_back(error(pi, pi2, i));	     		//progressive error of blocks
  } 
  
 Printfile_2( pi, err_pi, "Data3");
 
return 0;
}
