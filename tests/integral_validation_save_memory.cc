#include <iostream>
#include <string>
#include <vector>
#include <tuple>

using namespace std;

typedef struct{
  double m3pi;
  vector<double> diag;
  vector<string> names;
  vector<tuple<uint, uint, double, double> > all;
  pair<double, double> special;
} all_data;

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

void PrintVector(vector<string> v){
  for(uint i=0; i<v.size(); i++){
    cout << v[i] << endl;
  }
}

vector<all_data> ReadFiles(string fileprefix, string filesuffix, uint iwave1 = 0, uint iwave2 = 0, uint istart = 500, uint iend = 2490, uint ibin = 10){
  
  bool special = iwave1 != iwave2;
  if(istart<500 || istart>2490) istart = 500;
  if(iend>2490 || iend <istart) iend = 2490;
  if(ibin%10 != 0){
    cerr << "ibin is not dividable by 10. Corresponding files don't exist! Using ibin=10 instead" << endl;
    ibin = 10;
  }
  vector<all_data> output;
  for(uint ii=istart; ii<=iend; ii+=ibin){
    string filename = fileprefix + fourdigits(ii) + "_" + fourdigits(ii+10) + filesuffix;
    // cout << files[files.size()-1] << endl;
    all_data data;
    data.diag.resize(260);
    FILE* file = fopen(filename.c_str(),"r");
    char s1[20], s2[20], s3[20], s4[20];
    while(true){//find start of wave names
      fscanf(file, "%s", s1);
      if((string) s1 == "'FLAT") break;
    }
    data.names.push_back("FLAT");
    fscanf(file, "%s", s1);
    while((string) s1 == "'"){
      fscanf(file, "%s %s %s %s", s1, s2, s3, s4);
      if((string) s1 == "//") break;
      data.names.push_back((string) s1 + s2 + s3 + s4);
      fscanf(file, "%s", s1);
    }
    // cout << s1 << (((string)s1 == "//")?"true":"false") << endl;
    while(true){ //find starting point of the matrix
      fscanf(file, "%s", s1);
      // cout << 1 << endl;
      if((string)s1 == "normalisation") break; //second last word before the matrix starts
    };
    fscanf(file, "%s", s1); //still need to read the word "integrals"
    // cout << s1 << endl;
    for(uint i=0; i<260; i++){
      for(uint j=0; j<260; j++){
        double re, im;
        fscanf(file, " ( %lf, %lf)", &re, &im);
        // cout << re << ", " << im << endl;
        if(i==j){
          if(im != 0) cout << "ERROR: diagonal element is not ZERO!!!" << endl;
          data.diag[i] = re;
        }
        if(special){
          if(i == iwave1-1 && j == iwave2-1){
            data.special = make_pair(re, im);
            // cout << re << " " << im << endl;
          }
        }
        if(re==0 && im == 0) continue;
        data.all.push_back(make_tuple(i, j, re, im));
        // PrintTuple(data.all[data.all.size()-1]);
      }
    }
    
    // PrintVector(diag);
    
    // cout << "\%: " << values.size() << "/" << 260*260 << " = " << values.size()/(260.*260.) << endl;
    // double memory_needed = (2.*sizeof(uint) + 2.*sizeof(double))*values.size();
    // double memory_needed_before = 2.*sizeof(double)*260*260;
    // cout << "memory needed: " << memory_needed/1000000. << " MB\tmemory saved: " << (memory_needed_before - memory_needed)/1000000. << " MB" << endl;
    // cout << "\%memory: " << (double) memory_needed/(memory_needed_before - memory_needed) << endl;
    fclose(file);
    data.m3pi = ii/1000.;
    output.push_back(data);
    cout << "Reading: " << (ii-istart)*100/(iend-istart) << "\%             "<< "\r" << flush;
  }
  cout << "Reading finished.            " << endl;
  // PrintVector(output[0].diag);
  
  return output;
}


int main(){
  
  string fileprefix = "/mnt/data/compass/2008/integrals_Florian/PWANormIntegralsNAcc_";
  string filesuffix = "_0100_0113.txt";
  uint iwave1 = 13, iwave2 = 17;
  
  
  vector<all_data> output = ReadFiles(fileprefix, filesuffix, iwave1, iwave2);
  // PrintVector(output[0].names);
  // cout << output[0].names.size() << endl;
  string filename = "diag.txt";
  FILE* outputfile = fopen(filename.c_str(), "w");
  cout << "Writing to file..." << endl;
  fprintf(outputfile, "m3pi ");
  for(uint i=0; i<output.size(); i++){
    fprintf(outputfile, "%lf ", output[i].m3pi);
  }
  fprintf(outputfile, "\n\n");
  
  for(uint ii=0; ii<260; ii++){
    fprintf(outputfile, "%s ", output[0].names[ii].c_str());
    for(uint i=0; i<output.size(); i++){
      fprintf(outputfile, "%lf ", output[i].diag[ii]);
    }
    fprintf(outputfile, "\n");
  }
  fclose(outputfile);
  cout << "Writing to file '" << filename << "' completed successfully." << endl;
  
  if(iwave1 != iwave2){
    outputfile = fopen(("offdiag_" + to_string(iwave1) + "_" + to_string(iwave2) + ".txt").c_str(), "w");
    for(uint i=0; i<output.size(); i++){
      fprintf(outputfile, "%lf %lf %lf\n", output[i].m3pi, output[i].special.first, output[i].special.second);
    }
  }
  cout << "Writing to file 'offdiag_" << iwave1 << "_" << iwave2 << ".txt' completed successfully." << endl;
  
  // for(uint i=0; i< 260; i++){
  //   cout << i+1 << " " << output[0].names[i] << endl;
  // }
  return 0;
}
