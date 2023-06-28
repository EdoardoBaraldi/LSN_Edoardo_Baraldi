#include "mylib.h"

using namespace std;

const double a=1.0;   		//define lattice constant

const double n_step=100; 	//define nuber of steps

int main(){

 //cubic lattice
 vector<double> sum_prog(N, 0.0);  	// progressive values of r^2 of each block in lattice
 vector<double> sum2_prog(N, 0.0);	// progressive values of (r^2)^2 of each block in lattice
 vector<double> err_prog(N, 0.0);	// progressive values of r^2 values in lattice
 vector<double> sum_prog_cont(N, 0.0);  // progressive values of r^2 of each block continous
 vector<double> sum2_prog_cont(N, 0.0);	// progressive values of (r^2)^2 of each block continous
 vector<double> err_prog_cont(N, 0.0);	// progressive values of r^2 values continous
 
 
 
 vector<double> r2_avg(n_step+1, 0.0);  // values of sqrt(<r^2>) for each step in lattice	
 vector<double> r2_err(n_step+1,0.0);	// values of sqrt(<r^2>)'s error for each step in lattice
 vector<double> r2_avg_cont(n_step+1, 0.0);  	// values of sqrt(<r^2>) for each step continous	
 vector<double> r2_err_cont(n_step+1,0.0);	// values of sqrt(<r^2>)'s error for each step continous
 
 
 double sum=0.0; 			// sum of L values in each block in lattice
 double sum_cont=0.0;			// sum of L values in each block continous
 
 vector<vector<vector<double>>> traj_vector(N, vector<vector<double>>(L,vector<double>(3,0.0)));//vector of M simulations of random walks in lattice
 vector<vector<vector<double>>> traj_vector_cont(N, vector<vector<double>>(L,vector<double>(3,0.0)));//vector of M simulations of random walks continous
 
 Random Rand;			//random generator
 pray_the_sun(Rand); 		//initialize random generator
 
 for (int i = 1; i <= n_step; i++) {			//cycle of steps
  for (int block=0; block < N; block++) {		//cycle of N blocks
    for(int j=0; j<L; j++){				//cycle on the L elements of each block
     trajectory(Rand, traj_vector[block][j], a);	//random walk of M simulations for the lattice 
     trajectory_cont(Rand, traj_vector_cont[block][j], a);	//random walk of M simulations with continous position
     sum += pow(distance(traj_vector[block][j]),2);	//sum of r^2 of each random walk at each step in lattice
     sum_cont += pow(distance(traj_vector_cont[block][j]),2);  //sum of r^2 of each random walk at each step continous
    }
    if(block==0){
     sum_prog[block]=sum/double(L);				//mean of block 0 in lattice
     sum2_prog[block]=pow(sum/double(L),2);			//square mean of block 0 in lattice
     
     sum_prog_cont[block]=sum_cont/double(L);			//mean of block 0 continous
     sum2_prog_cont[block]=pow(sum_cont/double(L),2);		//square mean of block 0 continous
    }
    else{
     sum_prog[block]=sum_prog[block-1] + sum/double(L);  	  //progressive sum of r^2 for each block in lattice
     sum2_prog[block]=sum2_prog[block-1] + pow(sum/double(L),2);  //progressive sum of (r^2 )^2 for each block in lattice
     
     sum_prog_cont[block]=sum_prog_cont[block-1] + sum_cont/double(L);  	 //progressive sum of r^2 for each block continous
     sum2_prog_cont[block]=sum2_prog_cont[block-1] + pow(sum_cont/double(L),2);  //progressive sum of (r^2 )^2 for each block continous
    }
    
    sum=0;         					//reset sum to 0 for new block in lattice
    sum_cont=0;         				//reset sum to 0 for new block continous	
  }
  for(int block=0; block<N; block++){
   sum_prog[block] /= (block+1);				//progressive mean of blocks in lattice
   sum2_prog[block] /= (block+1);				//progressive square mean of blocks in lattice
   err_prog[block] = error (sum_prog, sum2_prog, block);	//progressive error of blocks in lattice
   
   sum_prog_cont[block] /= (block+1);					//progressive mean of blocks continous
   sum2_prog_cont[block] /= (block+1);					//progressive square mean of blocks continous
   err_prog_cont[block] = error (sum_prog_cont, sum2_prog_cont, block);	//progressive error of blocks continous
  }								
  
  r2_avg[i] = sqrt(sum_prog[N-1]);					//for each step i take the N-th progressive values in lattice
  r2_err[i] = sqrt(err_prog[N-1]) / (2*sqrt(sum_prog[N-1]));	//for each step i take the error propagation of the N-th progressive <r^2> values in lattice
  
  r2_avg_cont[i] = sqrt(sum_prog_cont[N-1]);				//for each step i take the N-th progressive values continous
  r2_err_cont[i] = sqrt(err_prog_cont[N-1]) / (2*sqrt(sum_prog_cont[N-1]));	//for each step i take the error propagation of the N-th progressive <r^2> values continous
  
  
  // i calculated the progressive values of r^2, but i want the sqrt of the mean progressive value therefore i have to propagate the error.
  
  fill(sum_prog.begin(), sum_prog.end(), 0.);		//fill the r^2 progressive vector with 0 in lattice
  fill(sum2_prog.begin(), sum2_prog.end(), 0.); 	//fill the (r^2)^2 progressive vector with 0 in lattice
  fill(err_prog.begin(), err_prog.end(), 0.);		//fill error progressive vector with 0 in lattice
  
  fill(sum_prog_cont.begin(), sum_prog_cont.end(), 0.);		//fill the r^2 progressive vector with 0 continous
  fill(sum2_prog_cont.begin(), sum2_prog_cont.end(), 0.); 	//fill the (r^2)^2 progressive vector with 0 continous
  fill(err_prog_cont.begin(), err_prog_cont.end(), 0.);		//fill error progressive vector with 0 continous
  
 }
 
 Printfile_2( r2_avg, r2_err, "Data2_1");		//print in a file the mean of each step in lattice with its own error
 Printfile_2( r2_avg_cont, r2_err_cont, "Data2_2");		//print in a file the mean of each step continous with its own error
    

return 0;
}
