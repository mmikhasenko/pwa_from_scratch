#include <iostream>
#include <string>
#include <vector>
#include <tuple>

using namespace std;


string fourdigits(uint i){
  if(i<1000){
    return "0" + to_string(i);
  }else{
    return to_string(i);
  }
}

void PrintTuple(tuple<uint, uint, double, double> t){
  cout << "(" << get<0>(t) << ", " << get<1>(t) << ", " << get<2>(t) << ", " << get<3>(t) << ")" << endl;
}

void PrintVector(vector<double> v){
  for(uint i=0; i<v.size(); i++){
    cout << v[i] << endl;
  }
}

int main(){
  
  string fileprefix = "/mnt/data/compass/2008/integrals_Florian/PWANormIntegralsNAcc_";
  string filesuffix = "_0100_0113.txt";
  
  vector<string> files;
  for(uint i=500; i<=2490; i+=10){
    files.push_back(fileprefix + fourdigits(i) + "_" + fourdigits(i+10) + filesuffix);
    // cout << files[files.size()-1] << endl;
  }
  vector<tuple<uint, uint, double, double> > values;
  vector<double> diag;
  diag.resize(260);
  FILE* file = fopen(files[0].c_str(),"r");
  char s1[20];
  fscanf(file, "%s", s1);
  // cout << s1 << (((string)s1 == "//")?"true":"false") << endl;
  while(true){
    fscanf(file, "%s", s1);
    // cout << 1 << endl;
    if((string)s1 == "normalisation") break;
  };
  fscanf(file, "%s", s1);
  // cout << s1 << endl;
  for(uint i=0; i<260; i++){
    for(uint j=0; j<260; j++){
      double re, im;
      fscanf(file, " ( %lf, %lf)", &re, &im);
      // cout << re << ", " << im << endl;
      if(i==j){
        if(im != 0) cout << "ERROR: diagonal element is not ZERO!!!" << endl;
        diag[i] = re;
      }
      if(re==0 && im == 0) continue;
      values.push_back(make_tuple(i, j, re, im));
      // PrintTuple(values[values.size()-1]);
    }
  }
  
  PrintVector(diag);
  
  // cout << "\%: " << values.size() << "/" << 260*260 << " = " << values.size()/(260.*260.) << endl;
  // double memory_needed = (2.*sizeof(uint) + 2.*sizeof(double))*values.size();
  // double memory_needed_before = 2.*sizeof(double)*260*260;
  // cout << "memory needed: " << memory_needed/1000000. << " MB\tmemory saved: " << (memory_needed_before - memory_needed)/1000000. << " MB" << endl;
  // cout << "\%memory: " << (double) memory_needed/(memory_needed_before - memory_needed) << endl;
}
