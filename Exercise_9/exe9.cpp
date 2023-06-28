#include "exe9.h"


//functions
void Initialize();
void Set_Positions_circle();
void Set_Positions_square();
void L1(vector<double> &cromo);
void L2(vector<double> &cromo);
void Reset_L();
double error_standard(double sum, double sum2, int N);
	

//main
int main() {
 Initialize();
 //Use L_1 
 Darwin.Set_Genome(n_cromosomes, n_cities, L1);
 for(int igen=1; igen<=generations; igen ++) {
  Darwin.Mutation_Permutation(Rand.Rannyu(), L1);
  Darwin.Mutation_Shift(Rand.Rannyu(), L1);
  Darwin.Mutation_Shuffle(Rand.Rannyu(), L1);
  Darwin.Mutation_Inversion(Rand.Rannyu(), L1);
  Darwin.Crossover(Rand.Rannyu(), L1);
  Reset_L();
  for(int irighe=0; irighe<int(n_cromosomes/2); irighe ++) {
   sum_L1 += Darwin.Get_L(irighe);
   sum_2_L1 += pow(Darwin.Get_L(irighe),2);
  }
  ofstream HalfL1;
  HalfL1.open("Results/averege_square_half_L1.data", ios::app);
  HalfL1 << igen << "," <<  Darwin.Get_L(0) << "," << sum_L1/(double(n_cromosomes)/2.) << "," << error_standard(sum_L1, sum_2_L1, int(n_cromosomes/2)) << endl;
  HalfL1.close();
 }
 ofstream BestL1;
 BestL1.open("Results/best_square_L1.data", ios::app);
 for(int ibest=0; ibest<n_cities; ibest ++) {
  BestL1 << Darwin.Get_L(0) << "," << Darwin.Get_City(0,ibest+1) << "," << Pos_city[int(Darwin.Get_City(0,ibest+1))-1][0] << "," << Pos_city[int(Darwin.Get_City(0,ibest+1))-1][1] << endl;
 }
 BestL1 << Darwin.Get_L(0) << "," << Darwin.Get_City(0,1) << "," << Pos_city[int(Darwin.Get_City(0,1))-1][0] << "," << Pos_city[int(Darwin.Get_City(0,1))-1][1] << endl;
 BestL1.close();

 //Use L_2
 /*Darwin.Set_Genome(n_cromosomes, n_cities, L2);
 for(int igen=1; igen<=generations; igen ++) {
  Darwin.Mutation_Permutation(Rand.Rannyu(), L2);
  Darwin.Mutation_Shift(Rand.Rannyu(), L2);
  Darwin.Mutation_Shuffle(Rand.Rannyu(), L2);
  Darwin.Mutation_Inversion(Rand.Rannyu(), L2);
  Darwin.Crossover(Rand.Rannyu(), L2);
  Reset_L();
  for(int irighe=0; irighe<int(n_cromosomes/2); irighe ++) {
   sum_L2 += Darwin.Get_L(irighe);
   sum_2_L2 += pow(Darwin.Get_L(irighe),2);
  }
  ofstream HalfL2;
  HalfL2.open("Results/averege_circ_half_L2.data", ios::app);
  HalfL2 << igen << "," <<  Darwin.Get_L(0) << "," << sum_L2/(double(n_cromosomes)/2.) << "," << error_standard(sum_L2, sum_2_L2, int(n_cromosomes/2)) << endl;
  HalfL2.close();
 }
 ofstream BestL2;
 BestL2.open("Results/best_circ_L2.data", ios::app);
 for(int ibest=0; ibest<n_cities; ibest ++) {
  BestL2 << Darwin.Get_L(0) << "," << Darwin.Get_City(0,ibest+1) << "," << Pos_city[int(Darwin.Get_City(0,ibest+1))-1][0] << "," << Pos_city[int(Darwin.Get_City(0,ibest+1))-1][1] << endl;
 }
 BestL2 << Darwin.Get_L(0) << "," << Darwin.Get_City(0,1) << "," << Pos_city[int(Darwin.Get_City(0,1))-1][0] << "," << Pos_city[int(Darwin.Get_City(0,1))-1][1] << endl;
 BestL2.close();
 */
 return 0;
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
   
   ReadInput >> generations;
   
   ReadInput >> n_cities;
   
   ReadInput >> n_cromosomes;

   ReadInput >> restart;
   
   if(restart) Seed.open("seed.out");
   else Seed.open("seed.in");
   Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   Rand.SetRandom(seed,p1,p2);
   Seed.close();
   
   ReadInput >> prob_mut[0] >> prob_mut[1] >> prob_mut[2] >> prob_mut[3];
   
   ReadInput >> prob_cross;
   
   ReadInput >> prob_selec;
   
   //set the parameters of Eugenetica class
   Darwin.Set_Probability(prob_cross, prob_mut, prob_selec);
   
   ReadInput >> Circ_Square;
   
   if(Circ_Square) Set_Positions_circle();
   else Set_Positions_square();
  
}

void Set_Positions_circle() {
   double theta;
   for(int i=0; i<n_cities; i++) {
    double r= 1.;
    theta = Rand.Rannyu(0., 2*M_PI);
    vector<double> pos = {r*cos(theta), r*sin(theta)};
    Pos_city.push_back(pos);
   }
}

void Set_Positions_square() {
   for(int i=0; i<n_cities; i++) {
    vector<double> pos = {Rand.Rannyu(0., 2), Rand.Rannyu(0., 2)};
    Pos_city.push_back(pos);
   }
}

void L1(vector<double> &cromo) {
    if(int(cromo.size()) == n_cities +1) {
     cromo[0] = 0.;
     double dist_x, dist_y;
     for(int i=1; i<(int(cromo.size())-1); i++) {
      dist_x = Pos_city[int(cromo[i])-1][0] - Pos_city[int(cromo[i+1])-1][0];
      dist_y = Pos_city[int(cromo[i])-1][1] - Pos_city[int(cromo[i+1])-1][1];
      cromo[0] += sqrt(pow(dist_x,2) + pow(dist_y,2));
     }
     dist_x = Pos_city[int(cromo[cromo.size()-1])-1][0] - Pos_city[int(cromo[1])-1][0];
     dist_y = Pos_city[int(cromo[cromo.size()-1])-1][1] - Pos_city[int(cromo[1])-1][1];
     cromo[0] += sqrt(pow(dist_x,2) + pow(dist_y,2));
    }
}

void L2(vector<double> &cromo) {
    if(int(cromo.size()) == n_cities +1) {
     cromo[0] = 0.;
     double dist_x, dist_y;
     for(int i=1; i<(int(cromo.size())-1); i++) {
      dist_x = Pos_city[int(cromo[i])-1][0] - Pos_city[int(cromo[i+1])-1][0];
      dist_y = Pos_city[int(cromo[i])-1][1] - Pos_city[int(cromo[i+1])-1][1];
      cromo[0] += (pow(dist_x,2) + pow(dist_y,2));
     }
     dist_x = Pos_city[int(cromo[cromo.size()-1])-1][0] - Pos_city[int(cromo[1])-1][0];
     dist_y = Pos_city[int(cromo[cromo.size()-1])-1][1] - Pos_city[int(cromo[1])-1][1];
     cromo[0] += (pow(dist_x,2) + pow(dist_y,2));
    }
}

void Reset_L() {
     sum_L1=0; 
     sum_2_L1=0;
     sum_L2=0;
     sum_2_L2=0;
     error_L1=0;
     error_L2=0;
}

double error_standard(double sum, double sum2, int N) {
     double err= sqrt(sum2/(double)N - pow(sum/(double)N,2));
     return err;
}









