#ifndef __ISING__
#define __ISING__

//libraries
#include "random.h"
#include "mylib.h"

using namespace std;

//simulations
Random Rand;		//Random generator
int seed[4];
int nstep, nblk;
double delta;		// size of the uniform proposal distribution 
int equilibration, restart, find_width; 
int steps_acc = 0;		//steps used to try the true width with acceptance rate of 50%

//averages
double blk_norm, blk_av, stima_H;
double sum, sum2, err;
double ACC=0;		//comulative acceptance rate
double accepted,attempted;

//Simluated annealing
int SA;
double beta; //initial temperature and beta_steps
double J, I; //metropolis steps to calculate H for each configuration
vector<double> param_old(4), param_new(4);
double x_new, x_old=0;  //position
double H, p, test;
double mu_old, mu_new, sigma_old, sigma_new, mu_initial, sigma_initial, H_old, H_err_old;

//hysto
int bin;
const int N_bins=1000;
double err_psi[N_bins], stima_psi[N_bins], glob_av_hysto[N_bins], glob_av2_hysto[N_bins], block_av_hysto[N_bins], walker_hysto[N_bins];
double width_bin=6./(double)N_bins;

//metropolis parameters
double a=0;	// acceptance probability 
double r;      // random number for acceptance

#endif

