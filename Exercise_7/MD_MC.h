/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid__
#define __fluid__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, iv, ik, it, ie, ip, iw;
double vtail, ptail, bin_size, nbins, sd;
double walker[m_props];
const double ref_temp = 1.2; // set the reference temperature 

// averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, stima_press, stima_kin, stima_etot, stima_temp;
double err_pot, err_press, err_kin, err_etot, err_temp;

//configuration
const int m_part=108;
double x[m_part],    y[m_part],    z[m_part];
double xold[m_part], yold[m_part], zold[m_part];
double vx[m_part],  vy[m_part],   vz[m_part];

// thermodynamical state
int npart;
double beta,temp,energy,vol,rho,box,rcut;
double ene_corr, press_corr; //tail corrections

// simulation
int iNVET, nstep, nblk, restart;
double delta;

//pigreco
const double pi=3.1415927;

//g(r)
int bin;
const int N_bins=500;
double delta_V[N_bins], err_g[N_bins], stima_g[N_bins];
double width_bin;

//equilibration
int index = 1;
double new_temp_eq;
bool U_on_N = false;	//parameter to control correltaion
bool acceptance = false; //parameter to control acceptance
bool equilib = false;	//parameter to equilibrate
int skip = 1000 ;	//how many steps i have to skip to have states equilibrated
double ACC = 0.;	//acceptance/attempted value

//functions
void Input(void);
void Reset(int);
void Accumulate(int,int);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);
double Force(int, int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
