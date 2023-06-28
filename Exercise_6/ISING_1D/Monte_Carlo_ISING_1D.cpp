/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

double min(const double a, const double b) {
	return (a <= b) ? a : b;
}

int main()
{ 
  Input(); //Inizialization
  string cycle_temp;
  cout << "Do you want to cycle over temperatures [0.5 ; 2]? (y/n) \n (if you're equilibrating it's not recommended)" << endl;
  cin >> cycle_temp;
  if( cycle_temp =="y") {	//cycle over temperatures to see the propery of the system vary with temperature
   cycle = true;         
   temp = 2.;		// initialize the temperature at 2
   beta = 1./temp;
   do {
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
     Reset(iblk);   //Reset block averages
     for(int istep=1; istep <= nstep; ++istep)
     {
      Move(metro);
      Measure();
      Accumulate(iblk, istep); //Update block averages
     }
     Averages(iblk);   //Print results for current block
    }
    ConfFinal(); //Write final configuration
    temp -= 0.05;
    beta = 1./temp;
   } while(temp >= 0.495);
  }
  
  
  else {
   cycle = false;
   for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
   {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(iblk, istep); //Update block averages
    }
    Averages(iblk);   //Print results for current block
   }
   ConfFinal(); //Write final configuration
  }

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();
   
   string starting;
   cout << "Do you want to equlibrate? \n (y/n)" << "\n" << endl;
   cin >> starting;
   
   if(starting =="y") {
    equilibration = true;
   }
   
   else { 
    equilibration = false;
   }
   
   string whatseed;
   cout <<"Do you want the original seed? \n (y/n)" << "\n" << endl;
   cin >> whatseed;
   
   if(whatseed =="y") {
    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
   }
   
   else {
    ifstream newseed("seed.out");
    newseed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    newseed.close();
   }
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep; 
  
  if(equilibration) {
   nblk = 10000;
   nstep = 1;
  }

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
	energy_old = Boltzmann(s[o], o);
	s[o] *= -1;
	energy_new = Boltzmann(s[o], o);
	
	attempted ++;
	
	double a= min(1., exp(-beta*(energy_new - energy_old)));
	if (rnd.Rannyu() >= a) {
	  s[o] *= -1.;
	  accepted ++;
	}
    }
    else //Gibbs sampling
    {
	energy_down = Boltzmann(-1,o);
	energy_up = Boltzmann(1,o);
	p = 1./(1. + exp(-beta*(energy_up - energy_down)));
	
	if(rnd.Rannyu() < p) {
	 s[o] = -1.;
	}
	else {
	 s[o] = 1.;
	}
	
	accepted ++;
	attempted ++;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
  
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(int iblock, int istep) //Update block averages
{
 if(equilibration) {
  for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
 }
 
 else {
  if(iblock = 1) {
   if(istep >= 250) {
    for(int i=0; i<n_props; ++i)
    {
     blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
   }
  }
  else {
    for(int i=0; i<n_props; ++i)
    {
     blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
  }
 }
 
}


void Averages(int iblk) //Print results for current block
{
  
  stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
  glob_av[iu]  += stima_u;
  glob_av2[iu] += stima_u*stima_u;
  err_u=Error(glob_av[iu],glob_av2[iu],iblk);
  
  stima_c = (beta * beta *(blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm,2)))/(double)nspin; //Heat capacity
  glob_av[ic]  += stima_c;
  glob_av2[ic] += stima_c*stima_c;
  err_c=Error(glob_av[ic],glob_av2[ic],iblk);
  
  stima_x = (beta * blk_av[ix]/blk_norm )/(double)nspin; //Susceptibility
  glob_av[ix]  += stima_x;
  glob_av2[ix] += stima_x*stima_x;
  err_x=Error(glob_av[ix],glob_av2[ix],iblk);
  
  stima_m = blk_av[im]/blk_norm /(double)nspin; //Magnetization
  glob_av[im]  += stima_m;
  glob_av2[im] += stima_m*stima_m;
  err_m=Error(glob_av[im],glob_av2[im],iblk);
  
  const int wd=12;
  
  if(cycle && !equilibration) {
   if(iblk == nblk) {
    ofstream Temp;
    Temp.open("gibbs_cycle_temp_02.dat", ios::app);
    // temperature -- energy -- energy_error -- heat capacity -- heat_error -- magnetisation -- mag_error -- chi -- chi_error
    Temp << temp << "," <<  glob_av[iu]/(double)iblk << "," << err_u << "," << glob_av[ic]/(double)iblk << "," << err_c << "," << glob_av[im]/(double)iblk << "," << err_m << "," << glob_av[ix]/(double)iblk << "," << err_x << endl;
    Temp.close();
   }
  }
    
  
  if(!cycle) {
   if(equilibration) {
    ofstream Equilib;
    Equilib.open("gibbs_equilib_05.dat", ios::app);
    Equilib << setw(wd) << iblk << setw(wd) << stima_u<< endl;
    Equilib.close();
    //cout << "----------------------------" << endl << endl;
   }
   
   else { 
    ofstream Ene, Heat, Mag, Chi;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output.ene.0",ios::app);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    Heat.open("output.heat.0",ios::app);
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();
    
    Chi.open("output.chi.0",ios::app);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();
    
    Mag.open("output.mag.0",ios::app);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Chi.close();
   }
  }
    
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
