#ifndef __GENETIC__
#define __GENETIC__
//libraries
#include "random.h"
#include "genetic_alg.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <vector>
#include <functional>
#include "mpi.h"
//#include "mylib.h"

using namespace std;

//random parameters
Random Rand;
int seed[4];
int restart;

//genetic algorithm
int Circ_Square;				//0 for circumference and 1 for square
int n_cities;
int n_cromosomes, generations;		//number of cities, cromosomes and generations
vector<vector<double>> Pos_city;	//vector of positions in 2D (x,y) for each city
vector<vector<std::string>> Name_City;
double prob_mut[4]; 				//probability of mutantions
double prob_cross, prob_selec;			//probability of crossover and selection
Eugenetica Darwin;				//Genetic class

//avereages
double sum_L1, sum_2_L1, sum_L2, sum_2_L2, error_L1, error_L2;	//sum, sum squared and error of L1 and L2

#endif // __GENETIC__

