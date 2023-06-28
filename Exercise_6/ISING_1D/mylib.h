#pragma once

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <assert.h>
#include <vector>
#include <algorithm>

#include "random.h"

using namespace std;
#define _USE_MATH_DEFINES

int M=pow(10,6); 	//Total Number of throws
int N=100; 		//Total Number of blocks
int L=int(M/N); 	//Number of throws in each block

double error(const vector<double> ave1, const vector<double> ave2, const int n){
	if(n==0){
		return 0;}
	else{
		return sqrt((ave2[n]-pow(ave1[n],2))/n);}
}
	

template <typename T> vector<T> ReadData( const char* filename) {
	vector<T> v;
	ifstream fin(filename);
	if (!fin) { 
		cerr <<"Non esiste il file"<< filename <<endl;
		exit(-1);
	}
	else {
		while(!fin.eof()) {
			T a;
			fin >> a;
			v.push_back(a);
			// assert( !(fin.eof() )&& "Lettura file finita" );	
		}
	}
	fin.close();
	return v;
}

template <typename T> void Printfile(const vector<T> &v, const char* filename){
	ofstream outf(filename);
	for(T elem: v){outf<<elem<<"\n";}
	outf.close();
}

template <typename T> void Printfile_2(const vector<T> &v1,const vector<T> &v2, const char* filename){
	ofstream outf(filename);
	if(v1.size() != v2.size()){
		cerr <<"Error: vectors with different sizes\n";}
	else{
		for(int i=0; i<int(v1.size()); i++){
			outf<< v1[i] <<","<< v2[i] <<"\n";}
	}
	outf.close();
}

template <typename T> void Printfile_3(const vector<T> &v1,const vector<T> &v2, const vector<T> &v3, const char* filename){
	ofstream outf(filename);
	if((v1.size() != v2.size()) && (v1.size() != v3.size()) && (v3.size() != v2.size())){
		cerr <<"Error: vectors with different sizes\n";}
	else{
		for(int i=0; i<int(v1.size()); i++){
			outf<< v1[i] <<","<< v2[i] <<","<< v3[i] << "\n";}
	}
	outf.close();
}	

void pray_the_sun(Random &rnd){
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;	
}

template <typename T> void trajectory(Random &rand, vector<T> &traj, const double a){
  if(traj.size() != 3) {
    cerr << "Error: dimension of vector is not 3" << endl;
    return;
  }

  int t = rand.Rannyu(0, 3);
  traj[t] += (rand.Rannyu() >= 0.5) ? a : -a;
}

template <typename T> void trajectory_cont(Random &rand, vector<T> &traj, const double a){
   if((traj.size()==3)) {
    double theta=rand.Rannyu(0.,2*M_PI);
    double phi=rand.Rannyu(0.,M_PI);
    traj[0] += a*sin(phi)*cos(theta);
    traj[1] += a*sin(phi)*sin(theta);
    traj[2] += a*cos(phi);
   }
   else {
    cerr << "Error: dimension of vector is not 3"<< endl;
   } 
}

template <typename T> double distance(const vector<T> &dr) {
   if(dr.size()==3) {
    const T dx = dr[0];
    const T dy = dr[1];
    const T dz = dr[2];

    return sqrt(dx*dx + dy*dy + dz*dz);
   }	
   else {
    cerr << "Error: dimension of vector is not 3"<< endl;
    return 0;
   }
}



     
	
	
	
	
	
	
	
	
	
	
	
	
