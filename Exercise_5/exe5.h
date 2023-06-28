#ifndef __ISING__
#define __ISING__

//libraries
#include "random.h"
#include "mylib.h"

//simulations
Random Rand;		//Random generator
int M=pow(10,6); 	//Total Number of throws
int N=100; 		//Total Number of blocks
int L=int(M/N); 	//Number of throws in each block

//possibilities
double accepted,attempted;
string gaus_yn;
string want_equilib;
string find_width;
//equilibration
bool equilibration = false;  //equilibration parameter
bool acceptance_rate = false;

int steps = 0;		//steps used to try the true width with acceptance rate of 50%
double ACC=0;	//comulative acceptance rate

//parameters simulations
double width = 3.; // width of the uniform proposal distribution 

//averages
double blk_norm, blk_av, stima_r;
double sum, sum2, err;

//metropolis parameters
double r_com=0;  // comulative distance
double a=0;	// acceptance probability 
double r;      // random number for acceptance

//positions
vector<double> pos_prime(3);   //new positions
vector<double> pos = {0.0, 0.0, 0.0}; //starting point (with psi_2_1_0 use 0,0,1)


#endif


















































