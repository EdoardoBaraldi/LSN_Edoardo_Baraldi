#include "exe10.h"

//functions
void Initialize(int rango);
void Set_Positions_square();
void Load_Cities();
void L1(vector<double> &cromo);
void Reset_L();
double error_standard(double sum, double sum2, int N);	

int main(int argc, char* argv[]) {

 int size, rank;
 MPI_Init(&argc, &argv);
 MPI_Comm_size(MPI_COMM_WORLD, &size);
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 
 Initialize(rank);
 
 //not comunicating
 /*for(int igen=1; igen<=generations; igen ++) {
  for(int iteration=1; iteration<=50; iteration++) {
   Darwin.Mutation_Permutation(Rand.Rannyu(), L1);
   Darwin.Mutation_Shift(Rand.Rannyu(), L1);
   Darwin.Mutation_Shuffle(Rand.Rannyu(), L1);
   Darwin.Mutation_Inversion(Rand.Rannyu(), L1);
   Darwin.Crossover(Rand.Rannyu(), L1);
   Reset_L();
  }
  for(int irighe=0; irighe<int(n_cromosomes/2); irighe ++) {
   sum_L1 += Darwin.Get_L(irighe);
   sum_2_L1 += pow(Darwin.Get_L(irighe),2);
  }
  ofstream HalfL1;
  HalfL1.open("not_comunitating/cities_rank"+ to_string(rank) +"nc.data", ios::app);
  HalfL1 << igen << "," <<  Darwin.Get_L(0) << "," << sum_L1/(double(n_cromosomes)/2.) << "," << error_standard(sum_L1, sum_2_L1, int(n_cromosomes/2)) << endl;
  HalfL1.close();
 }
 ofstream BestL1;
 BestL1.open("not_comunitating/best_city_rank"+ to_string(rank) +"nc.data", ios::app);
 for(int ibest=0; ibest<n_cities; ibest ++) {
  BestL1 << Darwin.Get_L(0) << "," << Darwin.Get_City(0,ibest+1) << "," << Pos_city[int(Darwin.Get_City(0,ibest+1))-1][0] << "," << Pos_city[int(Darwin.Get_City(0,ibest+1))-1][1] << endl;
 }
 BestL1 << Darwin.Get_L(0) << "," << Darwin.Get_City(0,1) << "," << Pos_city[int(Darwin.Get_City(0,1))-1][0] << "," << Pos_city[int(Darwin.Get_City(0,1))-1][1] << endl;
 BestL1.close();
 */
 //comunicating
 for(int igen=1; igen<=generations; igen ++) {
  for(int iteration=1; iteration<=50; iteration++) {
   Darwin.Mutation_Permutation(Rand.Rannyu(), L1);
   Darwin.Mutation_Shift(Rand.Rannyu(), L1);
   Darwin.Mutation_Shuffle(Rand.Rannyu(), L1);
   Darwin.Mutation_Inversion(Rand.Rannyu(), L1);
   Darwin.Crossover(Rand.Rannyu(), L1);
   Reset_L();
  }
  for(int irighe=0; irighe<int(n_cromosomes/2); irighe ++) {
   sum_L1 += Darwin.Get_L(irighe);
   sum_2_L1 += pow(Darwin.Get_L(irighe),2);
  }
  ofstream HalfL1;
  HalfL1.open("less_communicating/cities_rank"+ to_string(rank) +"lc.data", ios::app);
  HalfL1 << igen << "," <<  Darwin.Get_L(0) << "," << sum_L1/(double(n_cromosomes)/2.) << "," << error_standard(sum_L1, sum_2_L1, int(n_cromosomes/2)) << endl;
  HalfL1.close();
  //send and recive the best cromosomes
  vector<vector<double>> gatheredChromosomes;
  if (igen % 10000 == 0) {
   // Collect the chromosome from each process
   MPI_Barrier(MPI_COMM_WORLD);
   vector<double> bestchromosome(Darwin.Get_CromosomeSize());
   vector<double> receivedChromosomes(Darwin.Get_CromosomeSize() * size);  // Size è il numero totale di processi
   MPI_Gather(Darwin.Get_Cromosome(0).data(), Darwin.Get_CromosomeSize(), MPI_DOUBLE,
             receivedChromosomes.data(), Darwin.Get_CromosomeSize(), MPI_DOUBLE,
             0, MPI_COMM_WORLD);
   if (rank == 0) {
    // Aggiungi il cromosoma migliore di ogni processo al vettore di vettori gatheredChromosomes
    for (int i = 0; i < size; i++) {
     vector<double> chromosome(Darwin.Get_CromosomeSize());
     for (int j = 0; j < Darwin.Get_CromosomeSize(); j++) {
        chromosome[j] = receivedChromosomes[i * Darwin.Get_CromosomeSize() + j];
     }
     gatheredChromosomes.push_back(chromosome);
    }
   }
   // Only root process shuffles and distributes the chromosomes
   if (rank == 0) {
    random_shuffle(gatheredChromosomes.begin(), gatheredChromosomes.end());
   // Broadcast the shuffled chromosomes to all processes
    for (int i = 0; i < int(gatheredChromosomes.size()); i++) {
     MPI_Send(gatheredChromosomes[i].data(), Darwin.Get_CromosomeSize(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }
   }
   vector<double> receivedChromosome(Darwin.Get_CromosomeSize());
   MPI_Recv(receivedChromosome.data(), Darwin.Get_CromosomeSize(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   // Load the shuffled chromosomes in each process
   Darwin.Load_Cromosome(0, receivedChromosome);

  }
  // Clear the gatheredChromosomes vector
  gatheredChromosomes.clear();
 }
 ofstream BestL1;
 BestL1.open("less_communicating/best_city_rank"+ to_string(rank) +"lc.data", ios::app);
 for(int ibest=0; ibest<n_cities; ibest ++) {
  BestL1 << Darwin.Get_L(0) << "," << Darwin.Get_City(0,ibest+1) << "," << Pos_city[int(Darwin.Get_City(0,ibest+1))-1][0] << "," << Pos_city[int(Darwin.Get_City(0,ibest+1))-1][1] << endl;
 }
 BestL1 << Darwin.Get_L(0) << "," << Darwin.Get_City(0,1) << "," << Pos_city[int(Darwin.Get_City(0,1))-1][0] << "," << Pos_city[int(Darwin.Get_City(0,1))-1][1] << endl;
 BestL1.close();
 
 
 MPI_Finalize();
 return 0;
 
}


void Initialize(int rango) {

   ifstream ReadInput, Primes, Seed;
   Primes.open("Primes");
   int rowToRead = rango;
   int currentRow = 0;   	// Variable to keep track of the current row
   // Skip rows until the desired row is reached
   while(currentRow < rowToRead) {
    std::string line;
    std::getline(Primes, line);
    ++currentRow;
   }
   //Read seed for random numbers
   int p1, p2;
   Primes >> p1 >> p2 ;
   Primes.close();
    
   Seed.open("seed.in");
   Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   Rand.SetRandom(seed,p1,p2);
   Seed.close();
   
   //Read input informations
   ReadInput.open("input.in");
   
   ReadInput >> generations;
      
   ReadInput >> n_cromosomes;
   
   ReadInput >> prob_mut[0] >> prob_mut[1] >> prob_mut[2] >> prob_mut[3];
   
   ReadInput >> prob_cross;
   
   ReadInput >> prob_selec;
   
   //Set_Positions_square();
   
   Load_Cities();
   
   n_cities=int(Pos_city.size());
   
   //set the parameters of Eugenetica class
   Darwin.Set_Probability(prob_cross, prob_mut, prob_selec);
   Darwin.Set_Genome(n_cromosomes, n_cities, L1, rango);
   
   
}

void Load_Cities() {
   ifstream City;
   City.open("American_capitals.dat");
   double lon, lat;
   std::string line;
   std::string state;
   std::string name;
   bool isFirstRow = true;
   while (getline(City, line)) {
    istringstream Iss(line);
    if (isFirstRow) {
     isFirstRow = false;
     continue;
    }
    getline(Iss >> std::ws, state, '\t');
    getline(Iss >> std::ws, name, '\t');
    Iss >> lon >> lat; 
    Name_City.push_back({state,name});
    Pos_city.push_back({lon,lat});
   }
   City.close();
   ofstream Test_C;
   ofstream Test;
   Test.open("Test_Cities.dat");
   Test_C.open("Name_Cities.dat");
   for(int ibest=0; ibest<int(Pos_city.size()); ibest ++) {
    Test << ibest+1 << "," << Pos_city[ibest][0] << "," << Pos_city[ibest][1] << endl;
    Test_C << ibest+1 << "," << Name_City[ibest][0] << "," << Name_City[ibest][1] << endl;
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


void Reset_L() {
     sum_L1=0; 
     sum_2_L1=0;
     sum_L2=0;
     sum_2_L2=0;
     error_L1=0;
     error_L2=0;
}

double error_standard(double sum, double sum2, int N) {
     return sqrt(sum2/(double)N - pow(sum/(double)N,2));
}

