#include "mylib.h"

using namespace std;

int main(){

 vector<int> S_N{1,2,10,100};
 int S_realiz = int(pow(10,4));
 vector<double> sum_exp(S_realiz,0.);
 vector<double> sum_ran(S_realiz,0.);
 vector<double> sum_lor(S_realiz,0.);

 Random Rand;			//random generator
 pray_the_sun(Rand); 		//initialize random generator

 for(unsigned int i=0; i<S_N.size(); i++) {
  for(int j=0; j<S_realiz; j++) {
   for(int k=0; k<S_N[i]; k++){
    sum_exp[j]+=Rand.Exp(1.);
    sum_ran[j]+=Rand.Rannyu();
    sum_lor[j]+=Rand.Lorentz(0.,1.);
   }
  sum_exp[j]/=S_N[i];
  sum_ran[j]/=S_N[i];
  sum_lor[j]/=S_N[i];
  }
 string filename="Data2_" + to_string(S_N[i]);
 const char* filename_cstr= filename.c_str();
 Printfile_3(sum_ran, sum_exp, sum_lor, filename_cstr);
 fill(sum_ran.begin(), sum_ran.end(), 0.);
 fill(sum_lor.begin(), sum_lor.end(), 0.);
 fill(sum_exp.begin(), sum_exp.end(), 0.);
 } 
   
return 0;
}
