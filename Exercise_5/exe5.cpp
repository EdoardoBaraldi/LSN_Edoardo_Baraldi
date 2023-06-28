#include "exe5.h"

using namespace std;

//functions
double phi_1_0_0 (const vector<double> &pos);
double phi_2_1_0 ( const vector<double> &pos);
double T_unifrom(const vector<double> &pos_1, const vector<double> &pos_2, const double width);
double T_gaus(const vector<double> &pos_1, const vector<double> &pos_2, const double width);
void Print_a(void);
void Accomulate(int j);
void Move();
void ConfFinal();
void Reset(int i);
void Initialize();
void Print_equilib(int i);
void Averages(int i);


//main
int main(){
 
 Initialize();
 
 //if i want to find the correct acceptance i do a cycle until the accptance is between 0.49 and 0.51
 if(acceptance_rate) {				
  do {
   steps ++;
   accepted = 0.;
   attempted = 0.;
   ACC = 0.;
   fill(pos_prime.begin(), pos_prime.end(), 0.);
   fill(pos.begin(), pos.end(), 0.);
   for(int i=1; i<=N; i++){			//cycle on N blocks
    for(int j=1; j<=L; j++){			//cycle on the L elements of each block
     Move();
     Accomulate(j);
    }
    Averages(i);
   }
   width += 0.1;
  } while(((ACC/(double)N) >= 0.51) || ((ACC/(double)N) <= 0.49));
 }
  
 //if i want to sample the two probability distributions 
 else{				                //if i had already found the right width
  for(int i=1; i<=N; i++){			//cycle on N blocks
   Reset(i);	
   for(int j=1; j<=L; j++){			//cycle on the L elements of each block
    Move();
    if(equilibration) {   			//if i want to equilibrate i print in a file only the first 20000 steps
     Print_equilib(i);
    }
    else{
     Accomulate(j);
    }
    if(j%100 == 0) {
     ConfFinal();
    }
   }				//accomulate and update blco average and print the final result into a file
   if(!equilibration) {
    Averages(i);			//block averages and print into a file
   }
  }
 }
 
 return 0;
	
}


//functions' implementation
double phi_1_0_0 (const vector<double> &pos) {
	double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    	return exp(-2*r) / M_PI;
}

double phi_2_1_0 ( const vector<double> &pos) {
	double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
	if (r != 0) {
	 return 1./64.*(2/M_PI)*r*r*exp(-r)*pow(pos[2]/r,2);
	}
	else {
	 return 0;
	}
}

double T_unifrom(const vector<double> &pos_1, const vector<double> &pos_2, const double width) {
    double dist = 0.0;
    for (int i = 0; i<3; i++) {
        dist += pow(pos_1[i] - pos_2[i], 2);
    }
    if(sqrt(dist) <= width / 2.0) {
        return 1.0 / pow(width, 3);
    } else {
        return 0.0;
    }
}

double T_gaus(const vector<double> &pos_1, const vector<double> &pos_2, const double width) {
    double dist=0.0;
    for (int i = 0; i<3; i++) {
    	dist += pow(pos_1[i] - pos_2[i], 2);
    }
    double p = 1./((width/2.)*sqrt(2.*M_PI)) * exp(-dist/(2.*pow(width/2.,2)));
    return p;
}

//accomulation of block average
void Accomulate(int j) {
   blk_av = blk_av + r_com;
   blk_norm = blk_norm + 1.0;
}

//move to a new x
void Move() {
   for(int k=0; k<3; k++){			//cycle on the elemts of the position vector
    if(gaus_yn == "u") {	
     pos_prime[k] = Rand.Rannyu(pos[k] - width/2., pos[k] + width/2.);     // generate a proposal from a uniform distribution
    }
    if(gaus_yn == "g") {
     pos_prime[k] = Rand.Gauss(pos[k], width/2.);
    }
   }
   a = min(1.0, phi_2_1_0(pos_prime) / phi_2_1_0(pos)); // calculate acceptance probability with  probability function, moreover T is simmetric for both gaussian and uniform distribution so i can't simplify it the expression of a
   if (phi_2_1_0(pos) == 0) {
    a = 1.0;
   }
   attempted ++;
   r = Rand.Rannyu() ; // generate a uniform random number
   if (r <= a) {
    pos = pos_prime; // accept the move else reject and position remain the same 
    accepted ++;
   }
   r_com = distance(pos); //update r comulative
}

//reset variables
void Reset(int i) {
   if(i == 1) {
    sum = 0.;
    sum2 = 0.;
   }
   blk_av = 0.;
   blk_norm = 0;
   r_com=0;
   a=0;
   
   fill(pos_prime.begin(), pos_prime.end(), 0.);
   fill(pos.begin(), pos.end(), 0.);
}

//initialize the variable and also give the option to search the right acceptance rate and which T you want to use
void Initialize() {
   pray_the_sun(Rand); 		//initialize random generator
 
   cout <<"Which ditribution do you wnat to use: Gaussian or Uniform? \n (g/u)" <<endl;
   cin >> gaus_yn; 
   
   cout <<"Do you want to find the right width? \n (y/n)" <<endl;
   cin >> find_width;
   
   if(find_width =="y") {
    N = pow(10,6);
    L = 1;
    acceptance_rate=true;
   }
 
   cout <<"Do you want to equilibrate? \n (y/n) \t (it's not raccomanded if you want to find the right width)" <<endl;
   cin >> want_equilib;
 
   if((want_equilib =="y")) {equilibration=true;}
 
   if(equilibration) {
    N = 20000;
    L = 1;
   }  
}

//print final equilibration results in file
void Print_equilib(int i) {
   ofstream Equilib;
   if(gaus_yn == "g") {Equilib.open("equilib_2_1_0_gaus.dat", ios::app);}
   else {Equilib.open("equilib_2_1_0_uniform.dat", ios::app);}
   Equilib << i <<","<<distance(pos)<<endl;
   Equilib.close();
}

//evaluete the final r with the two probability distributions and print their block average in a file
void Averages(int i) {
   stima_r = blk_av/blk_norm;
   sum += stima_r;				//progressive mean of blocks
   sum2 += stima_r*stima_r;			//progressive square mean of blocks
   err = Error(sum, sum2, i);			//progressive error of blocks 
   
   ACC += accepted/attempted;
   
   if(acceptance_rate) {
    ofstream AccFile;
    if(gaus_yn == "g") {AccFile.open("acceptance_folder/accept_2_1_0_gaus.dat", ios::app);}
    else {AccFile.open("acceptance_folder/accept_2_1_0_uniform.dat", ios::app);}
     //Acceptance
    if( i==N ) {
     AccFile << steps <<  "," << ACC/(double)i << "," << width << endl;
     cout << steps <<  "," << ACC/(double)i << "," << width << endl;
    }
    AccFile.close();
   }
   
   else {
    ofstream R;
    if(gaus_yn == "u") {
     R.open("Data_1_0_0_uniform.dat", ios::app);
     R << i << "," << stima_r << "," << sum/double(i) << "," << err << endl;
     R.close();
     //print in a file 
    } 
    if(gaus_yn == "g") {
     R.open("Data_1_0_0_gaus.dat", ios::app);
     R << i << "," << stima_r << "," << sum/double(i) << "," << err << endl;
     R.close();
    }
   }
}	

//print the final configuration of the points distribution
void ConfFinal(void) {
   if(gaus_yn == "u") {
    ofstream WriteConf;
    WriteConf.open("config_1_0_0.final", ios::app);
    WriteConf << pos[0] << "," << pos[1] << "," << pos[2] <<endl;
    WriteConf.close();
   }
}  












