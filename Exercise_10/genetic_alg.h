#ifndef __Eugenetica__
#define __Eugenetica__

#include <vector>
#include <functional>
#include "random.h"


class Eugenetica {

private:
	std::vector<std::vector<double>> Genome;
	double p_c; 	//probability crossover
	double p_m1, p_m2, p_m3, p_m4;	//probability mutation
	double p_s;	//probability selection
	Random Nature;	//Random Generator
	int columns, rows;
	
	
public:
	// constructors
  	Eugenetica();
  	// destructor
  	~Eugenetica();
	void Set_Probability(double, double*, double);
	void Set_Genome(int, int, std::function<void(std::vector<double>&)>, int );
	void Reset_Genome();
	void Evaluate_L(std::function<void(std::vector<double>&)>);
	void Sorting();
	int Selection();
	void Check(int);
	void Mutation_Permutation(double, std::function<void(std::vector<double>&)>);
	void Mutation_Shift(double, std::function<void(std::vector<double>&)>);
	void Mutation_Shuffle(double, std::function<void(std::vector<double>&)>);
	void Mutation_Inversion(double, std::function<void(std::vector<double>&)>);
	void Crossover(double, std::function<void(std::vector<double>&)>);
	double Get_L(int);
	double Get_City(int, int);
	std::vector<double> Get_Cromosome(int);
	void Load_Cromosome(int, std::vector<double>);
	int Get_CromosomeSize();

};
	
#endif // __Eugenetica__





