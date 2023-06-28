#include "exe8.h"

//functions
double phi_T (const double x, const double mu, const double sigma);
double T_unifrom(const double &pos_1, const double &pos_2, const double width);
double Boltzamann(const double old_en, const double new_en, const double temp);
double V(const double x);
double K(const double x, const double mu, const double sigma);
double E(const double x, const double mu, const double sigma);
void Initialize();
void Move();
void Accumulate();
void ConfFinal();
void Reset(int iblk);
void Averages(int iblk);
void Print_equilib(int iblk);
vector<double> Get_Metropolis(int L, double mu, double sigma );
void Print_beta(double ibeta);


//main
int main(){
 
 Initialize();
 
 //simulated annealing to find the right parameters 
 if(SA) {
  mu_old=mu_initial;
  sigma_old=sigma_initial;
  
  for(double ibeta=beta; ibeta<=pow(10,3); ibeta+=(2.5*beta)) {
   cout << 1/ibeta << endl;
   for(int istep_beta=(int(1./ibeta*1000)); istep_beta>0; istep_beta--) {
    param_old = Get_Metropolis( L, mu_old, sigma_old );
    H_old=param_old[0];
    H_err_old=param_old[1];
    
    mu_new = fabs(Rand.Rannyu(mu_old - 1./ibeta, mu_old + 1./ibeta));
    sigma_new = fabs(Rand.Rannyu(sigma_old - 1./ibeta, sigma_old + 1./ibeta));
    
    param_new = Get_Metropolis(L, mu_new, sigma_new ); 
    test = min(1,exp(-ibeta*(param_new[0] -param_old[0])));
    if(test >= Rand.Rannyu()) {
     mu_old = mu_new;
     sigma_old = sigma_new;
     H_old = param_new[0];
     H_err_old = param_new[1];
    }
   }
   Print_beta(ibeta);
  }
  
  //after i've optimized the mu and sigma parameters i have to determinate the right delta for an acceptance of 50%
  mu_initial=mu_old;
  sigma_initial=sigma_old;
  
  double right_nblk=nblk;	//i create two support parameters for number of blocks and steps
  double right_nstep=nstep;
  
  find_width=1;
  delta=1.0;
  do {
   steps_acc ++;
   accepted = 0.;
   attempted = 0.;
   ACC = 0.;
   x_old=0;
   for(int iblk=1; iblk<=nblk; iblk++){			//cycle on N blocks
    for(int istep=1; istep<=nstep; istep++){			//cycle on the L elements of each block
     Move();
     Accumulate();
    }
    Averages(iblk);
   }
   delta += 0.1;
  } while(((ACC/(double)nblk) >= 0.52) || ((ACC/(double)nblk) <= 0.48));
  
  //after found the right delta for the Metropolis Algorithm i can evaluete the Energy 
  nblk=right_nblk;
  nstep=right_nstep;
  
  find_width=0;		//turn off the parameter to find the right delta
  
  for(int iblk=1; iblk<=nblk; iblk++){			//cycle on N blocks
   Reset(iblk);	
   for(int istep=1; istep<=nstep; istep++){		//cycle on the L elements of each block
    Move();
    Accumulate();
    if(istep%100 == 0) {
     ConfFinal();
    }
   }				//accomulate and update block average and print the final result into a file
   Averages(iblk);			//block averages and print into a file
  }
 }
 
 //no simualation annealing for testing the code
 else {
  if(find_width) {
   do {
    steps_acc ++;
    accepted = 0.;
    attempted = 0.;
    ACC = 0.;
    x_old=0;
    for(int iblk=1; iblk<=nblk; iblk++){			//cycle on N blocks
     for(int istep=1; istep<=nstep; istep++){			//cycle on the L elements of each block
      Move();
      Accumulate();
     }
     Averages(iblk);
    }
    delta += 0.5;
   } while(((ACC/(double)nblk) >= 0.52) || ((ACC/(double)nblk) <= 0.48));
  }
   
  else{				                //if i had already found the right width
   for(int iblk=1; iblk<=nblk; iblk++){		//cycle on N blocks
    Reset(iblk);	
    for(int istep=1; istep<=nstep; istep++){	//cycle on the L elements of each block
     Move();
     Accumulate();
     if(istep%100 == 0) {
      ConfFinal();
     }
    }						//accomulate and update blco average and print the final result into a file
    if(!equilibration) {
     Averages(iblk);				//block averages and print into a file
    }
    if(equilibration) {   			//if i want to equilibrate i print in a file only the first 20000 steps
     Print_equilib(iblk);
    }
   }
  }
 }

 return 0;
}



//functions
double Boltzamann(const double old_en, const double new_en, const double temp) {
	return exp(-1./temp*(new_en - old_en));
}

double phi_T (const double x, const double mu, const double sigma) {
    	return (exp(-pow(x-mu,2)/(2*sigma*sigma)) + exp(-pow(x+mu,2)/(2*sigma*sigma)));
}


double T_unifrom(const double &pos_1, const double &pos_2, const double width) {
    double dist = 0.0;
    dist += pow(pos_1 - pos_2, 2);
    if(sqrt(dist) <= width / 2.0) {
        return 1.0 / pow(width, 3);
    } else {
        return 0.0;
    }
}

double V(const double x) {
     return ((x*x -2.5)*pow(x,2.));
}

double K(const double x, const double mu, const double sigma) {
      return (-1./(sigma*sigma)*(exp(-pow(x-mu,2)/(2*sigma*sigma))*(1-pow((x-mu)/sigma,2)) + exp(-pow(x+mu,2)/(2*sigma*sigma))*(1-pow((x+mu)/sigma,2)))/phi_T(x,mu,sigma));
}

double E(const double x, const double mu, const double sigma) {
      return V(x) - 0.5 * K(x,mu,sigma);
}


void Initialize() {
   ifstream ReadInput, Primes, Seed;
   
   //Read seed for random numbers
   int p1, p2;
   Primes.open("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

//Read input informations
   ReadInput.open("input.in");
   
   ReadInput >> nblk;

   ReadInput >> nstep;  

   ReadInput >>restart;
   
   if(restart) Seed.open("seed.out");
   else Seed.open("seed.in");
   Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   Rand.SetRandom(seed,p1,p2);
   Seed.close();
   
   ReadInput >> find_width;
   
   if(find_width) {
    nblk = pow(10,6);
    nstep = 1;
   }
 
   ReadInput >> equilibration;
 
   if(equilibration) {
    nblk = 20000;
    nstep = 1;
   }
   
   ReadInput >> SA;
   
   ReadInput >> beta;
   
   ReadInput >> mu_initial;
   
   ReadInput >> sigma_initial;
   
   ReadInput >> J;
   
   ReadInput >> I;
   
   ReadInput >> delta;
}


void Move() {
   	
   x_new =Rand.Rannyu(x_old - delta/2., x_old + delta/2.);     // generate a proposal from a uniform distribution
   
   a = min(1.0, pow(fabs(phi_T(x_new,mu_initial,sigma_initial)),2) / pow(fabs(phi_T(x_old,mu_initial,sigma_initial)),2)); // calculate acceptance probability with  probability function, moreover T is simmetric for both gaussian and uniform distribution so i can't simplify it the expression of a
   attempted ++;
   r = Rand.Rannyu() ; // generate a uniform random number
   if (r <= a) {
    x_old = x_new; // accept the move else reject and position remain the same 
    accepted ++;
   }
   H = E(x_old,mu_initial,sigma_initial); //update h value
   
   //hysto fill
   if((x_old <= 3.) && (x_old >=-3)) {
    walker_hysto[int((x_old+3.)/width_bin)]+=1;
   }
}

void Accumulate() //Update block averages
{
  blk_av = blk_av + H;
  blk_norm = blk_norm + 1.0; 
  
  for(int i=0; i<N_bins; i++) {
   block_av_hysto[i] += walker_hysto[i];
  } 
}

void ConfFinal(void) {
    ofstream WriteConf;
    WriteConf.open("config_1_0_0.final", ios::app);
    WriteConf << x_old <<endl;
    WriteConf.close();
} 

void Reset(int iblk) {
   if(iblk == 1) {
    sum = 0.;
    sum2 = 0.;
    for(int i=0; i<N_bins; i++) {
     glob_av_hysto[i]=0.;
     glob_av2_hysto[i]=0.;
    }
   }
   for(int i=0; i<N_bins; i++) { walker_hysto[i] = 0.;}
   blk_av = 0.;
   blk_norm = 0;
   for(int i=0; i<N_bins; i++) {
     block_av_hysto[i]=0.;
   }
   a=0;
   if(!equilibration) {
    x_new=0;
    x_old=0;
   }
}

void Averages(int iblk) {
   stima_H = blk_av/blk_norm;
   sum += stima_H;				//progressive mean of blocks
   sum2 += stima_H*stima_H;			//progressive square mean of blocks
   err = Error(sum, sum2, iblk);			//progressive error of blocks 
   
   for(int i=0; i < N_bins; i++) {
      stima_psi[i] = block_av_hysto[i]/blk_norm;
      glob_av_hysto[i] += stima_psi[i];
      glob_av2_hysto[i] += stima_psi[i]*stima_psi[i];
      err_psi[i] = Error(glob_av_hysto[i],glob_av2_hysto[i],iblk);
   }
   
   ACC += accepted/attempted;
   
   if(find_width) {
    ofstream AccFile;
    AccFile.open("acceptance_folder/accept_H_uniform.dat", ios::app);
     //Acceptance
    if( iblk==nblk ) {
     AccFile << steps_acc <<  "," << ACC/(double)iblk << "," << delta << "," << endl;
     cout << steps_acc <<  "," << ACC/(double)iblk << "," << delta << endl;
    }
    AccFile.close();
   }
   
   else {
    ofstream R, Hysto;
    if(SA) {
     R.open("results/H_uniform_optimized_test4.dat", ios::app);
     Hysto.open("results/psi_hysto_test4.dat", ios::app);
     if(iblk == nblk) {
      for (int i=0; i < N_bins; i++) {
       Hysto << (i*width_bin -3.) << "," << glob_av_hysto[i]/(double)iblk << "," << err_psi[i] << endl;
      }
     }
    }  
    else {R.open("results/H_uniform_test1.dat", ios::app);}
    R << iblk << "," << stima_H << "," << sum/(double)iblk << "," << err << endl;
    R.close();
    Hysto.close();
   }
}

void Print_equilib(int iblk) {
   stima_H = blk_av/blk_norm;
   sum += stima_H;				//progressive mean of blocks
   ofstream Equilib;
   Equilib.open("equilib_folder/equilib_H_uniform.dat", ios::app);
   Equilib << iblk <<","<< sum/(double)iblk <<endl;
   Equilib.close();
}

vector<double> Get_Metropolis(int L, double mu, double sigma ) {
    for(int i=1; i<=I; i++) {
     Reset(i);
     for(int j=1; j<=J; j++) {
      x_new =Rand.Rannyu(x_old - delta/2., x_old + delta/2.);     // generate a proposal from a uniform distribution
      double acc=min(1,pow(fabs(phi_T(x_new,mu,sigma)),2) / pow(fabs(phi_T(x_old,mu,sigma)),2));
      if(Rand.Rannyu() <= acc) {
      //Update
       x_old= x_new;
      }
      H = E(x_old,mu,sigma); //update h value
      
      //accumulate
      blk_av = blk_av + H;
      blk_norm = blk_norm + 1.0;
     }
     stima_H = blk_av/blk_norm;
     sum += stima_H;				//progressive mean of blocks
     sum2 += stima_H*stima_H;			//progressive square mean of blocks
     err = Error(sum, sum2, i);		//progressive error of blocks 
    }
    
    vector<double> get_information(5);
    
    get_information[0]=sum/(double)I;		//H value
    get_information[1]=err;			//H err
    get_information[2]=mu;		    	//mu value
    get_information[3]=sigma;			//sigma value
    
    return get_information;
}

void Print_beta(double ibeta) {
   ofstream Beta;
   Beta.open("simulated_annealing/annealing_H_uniform_test4.dat", ios::app);
   Beta << ibeta <<","<< H_old <<","<< H_err_old <<","<< mu_old <<","<< sigma_old <<endl;
   Beta.close();
}
