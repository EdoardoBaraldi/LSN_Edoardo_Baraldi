#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <vector>
#include <functional>
#include "genetic_alg.h"


using namespace std;

Eugenetica :: Eugenetica(){}

Eugenetica :: ~Eugenetica(){}

void Eugenetica :: Set_Probability(double cross, double* mut, double sel) {
	p_c=cross;
	p_m1=mut[0];
	p_m2=mut[1];
	p_m3=mut[2];
	p_m4=mut[3];
	p_s=sel;
	return;
}

void Eugenetica :: Set_Genome(int righe, int colonne, function<void(vector<double>&)> L, int rango) {
        //initialize Random generator
        ifstream Primes, Seed;
        //Read seed for random numbers
        int p1, p2;
        Primes.open("Primes");
        int rowToRead = rango;  	// Row number to read 
        int currentRow = 0;   	// Variable to keep track of the current row
        // Skip rows until the desired row is reached
        while (currentRow < rowToRead) {
          std::string line;
          std::getline(Primes, line);
          ++currentRow;
        }
        Primes >> p1 >> p2 ;
        Primes.close();
        Seed.open("seed.in");
        int seed[4];
        Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        Nature.SetRandom(seed,p1,p2);
        Seed.close();
	rows=righe;
	columns=colonne+1;
	vector<double> firstrow;
	for(int i=0; i<columns; i++) {
	 firstrow.push_back(i);
	}
	Genome.push_back(firstrow);
	for(int j=1; j<rows; j++) {
	 Genome.push_back(Genome[0]);
	 random_shuffle(Genome[j].begin()+2, Genome[j].end());
	}
	for(int j=0; j<rows; j++) {
	 L(Genome[j]);
	}
	return;
}

void Eugenetica :: Reset_Genome() {
	Genome.clear();
}	 

void Eugenetica :: Evaluate_L(function<void(vector<double>&)> L) {
	for(int j=0; j<rows; j++) {
	 L(Genome[j]);
	}
	return;
}

void Eugenetica :: Sorting() {
    	sort(Genome.begin(), Genome.end(), [](const vector<double>& row1, const vector<double>& row2) {return row1[0] < row2[0];});
    	return;
}

int Eugenetica :: Selection() {
	Sorting();
	return (int(rows*pow(Nature.Rannyu(),p_s)));	
}	

void Eugenetica :: Mutation_Permutation(double prob, function<void(vector<double>&)> L) {
	if(prob < p_m1) { 
	 int best=Selection();
	 Sorting();
	 vector<double> row_best=Genome[best];
	 int m = int(Nature.Rannyu(1.,columns-2));
	 random_shuffle(row_best.begin() + 1+m, row_best.begin() +1+m+2);
	 Genome[Genome.size()-1] = row_best;
	 Check(Genome.size()-1);
	 Evaluate_L(L);
	 Sorting();
	}
	return;
}

void Eugenetica :: Mutation_Shift(double prob, function<void(vector<double>&)> L) {
	if(prob < p_m2) { 
	 int best=Selection();
	 Sorting();
	 vector<double> row_best=Genome[best];
	 int m = int(Nature.Rannyu(1.,int((columns-2)/2) +1));
	 int n = int(Nature.Rannyu(double(m)+1,double(columns-2-m)+1));
	 vector<double> first;
	 for(int i=1; i<=m; i++) {
	  first.push_back(row_best[i+1]);
	 }
	 vector<double> last;
	 for(int i=1; i<=m; i++) {
	  last.push_back(row_best[i+n+1]);
	 }
	 for(int i=1; i<=m; i++) {
	  row_best[i+1]=last[i-1];
	  row_best[i+1+n]=first[i-1];
	 }
	 Genome[Genome.size()-1] = row_best;
	 Check(Genome.size()-1);
	 Evaluate_L(L);
	 Sorting();
	}
	return;
}

void Eugenetica :: Mutation_Shuffle(double prob, function<void(vector<double>&)> L) {
	if(prob < p_m3) {
	 int best=Selection();
	 Sorting();
	 vector<double> row_best=Genome[best];
	 int m = int(Nature.Rannyu(1.,double(columns)-1.));
	 int n = int(Nature.Rannyu(1, double(columns-m)-1.));
	 random_shuffle(row_best.begin() + 1+m, row_best.begin() + 1+m+n);
	 Genome[Genome.size()-1] = row_best;
	 Check(Genome.size()-1);
	 Evaluate_L(L);
	 Sorting();
	}
	return;
}
	
void Eugenetica :: Mutation_Inversion(double prob, function<void(vector<double>&)> L) {
	if(prob < p_m4) {
	 int m = int(Nature.Rannyu(1.,double(columns)-1.));
	 int best=Selection();
	 Sorting();
	 vector<double> row_best=Genome[best];
	 reverse(row_best.begin() + 2, row_best.begin() + 2+m);
	 Genome[Genome.size()-1] = row_best;
	 Check(Genome.size()-1);
	 Evaluate_L(L);
	 Sorting();
	}
	return;
}

void Eugenetica :: Check(int row) {
	vector<double> copy_genome(columns-1);
	copy(Genome[row].begin() + 1, Genome[row].end(), copy_genome.begin());
	sort(copy_genome.begin(), copy_genome.end());
	if((adjacent_find(copy_genome.begin(), copy_genome.end(), std::equal_to<double>()) !=copy_genome.end()) &&  (std::accumulate(copy_genome.begin(), copy_genome.end(), 0.0) != (columns*(1+columns)/2))) {
	 copy_genome.clear();
	 cerr <<"There is a problem in the code, check your code" << endl;
	 return;
	}
	else {
	 copy_genome.clear();
	 return;
	}
}

void Eugenetica :: Crossover(double prob, function<void(vector<double>&)> L) {
	if(prob < p_c) {
	 int m = int(Nature.Rannyu(1.,double(columns)-1));
	 int father=Selection();
	 int mother=Selection();
	 vector<double> row_father=Genome[father];
	 vector<double> row_mother=Genome[mother];
	 vector<double> row_son;
	 vector<double> row_daughter;
	 copy(row_father.begin(), row_father.begin() + m+2, back_inserter(row_son));
	 copy(row_mother.begin(), row_mother.begin() + m+2, back_inserter(row_daughter));
	 for(int i=2; i<=columns-1; i++) {
	  if(all_of(row_son.begin(), row_son.end(),[row_mother,i](double elem){return elem!=row_mother[i];})) {
	   row_son.push_back(row_mother[i]);
	  }
	  if(all_of(row_daughter.begin(), row_daughter.end(),[row_father,i](double elem){ return elem!=row_father[i];})) {
	   row_daughter.push_back(row_father[i]);
	  }
	 }
	 Genome[Genome.size()-1] = row_son;
	 Check(int(Genome.size()-1));
	 Genome[Genome.size()-2] = row_daughter;
	 Check(int(Genome.size()-2));
	 Evaluate_L(L);
	 Sorting();
	}
	return;
}
		 
double Eugenetica :: Get_L(int riga) {
        return Genome[riga][0];
}

double Eugenetica :: Get_City(int riga, int colonna) {
        return Genome[riga][colonna];
}	
	 
vector<double> Eugenetica :: Get_Cromosome(int riga) {
	return Genome[riga];
}	

void Eugenetica :: Load_Cromosome(int row, std::vector<double> Cromosoma) {
	Genome[row] = Cromosoma;
	Sorting();
	return;
}

int Eugenetica :: Get_CromosomeSize(){
	return int(Genome[0].size());
}
        

    
		

