#include <fstream>
#include <string>
using namespace std;
int main(){
  ifstream is("run_list.txt");
  ofstream os("hadd.sh");
  os << " declare -a arr_10 arr_30 arr_50 arr_80 arr_100 arr_150 arr_mu\n";
  string runN;
  
  while(true){
    is >> runN;
    int N = 0;;
    int max_digi = runN.length();
    for( int i = 0 ; i < max_digi ; ++i){
      N *= 10;
      N += runN[i] - int('0');}
    if( N == 21 || N == 22) os << "arr_80+=( HexaOutput_" << N << ".root )\n";
    if( N == 23 ) os << "arr_50+=( HexaOutput_" << N << ".root )\n";
    if( N >= 24 && N <= 87) os << "arr_100+=( HexaOutput_" << N << ".root )\n";
    if( N == 88) os << "arr_mu+=( HexaOutput_" << N << ".root )\n";
    if( N >= 91 && N <= 131) os << "arr_80+=( HexaOutput_" << N << ".root )\n";
    if( N >= 132 && N <= 161) os << "arr_mu+=( HexaOutput_"<< N << ".root )\n";
    if( N >= 163 && N <= 204){
      if(!(N >= 169 && N <= 176))
	os << "arr_50+=( HexaOutput_" << N << ".root )\n";}
    if( N >= 205 && N <= 225 || N == 265) os << "arr_mu+=( HexaOutput_"<< N << ".root )\n";
    if( N >= 235 && N <= 263) os << "arr_100+=( HexaOutput_" << N << ".root )\n";
    if( N >= 268 && N <= 323) os << "arr_30+=( HexaOutput_" << N << ".root )\n";
    if( N >= 326 && N <= 363) os << "arr_10+=( HexaOutput_" << N << ".root )\n";
    if( N >= 365 && N <= 371) os << "arr_80+=( HexaOutput_" << N << ".root )\n";
    if( N >= 372 && N <= 376) os << "arr_50+=( HexaOutput_" << N << ".root )\n";
    if( N == 378 ) os << "arr_100+=( HexaOutput_" << N << ".root )\n";
    if( N >= 380 && N <=381 ) os << "arr_150+=( HexaOutput_" << N << ".root )\n";
    if( N >= 387 && N <= 392) os << "arr_mu+=( HexaOutput_"<< N << ".root )\n";
    if( N >= 393 && N <= 405) os << "arr_100+=( HexaOutput_" << N << ".root )\n";
    if( N >= 407 && N <=418 ) os << "arr_150+=( HexaOutput_" << N << ".root )\n";


    if(is.eof()) break;
  }
  
  
  return 0;}
