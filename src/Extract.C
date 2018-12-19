#include "TDirectory.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include <functional>


double invmasssq(double p2[ ])
{
  return p2[3]*p2[3]-p2[2]*p2[2]-p2[1]*p2[1]-p2[0]*p2[0];
}

double energy(double p2[ ], double masssq)
{
  return sqrt(masssq + (p2[2]*p2[2]+p2[1]*p2[1]+p2[0]*p2[0]));
}

double costheta(double p[])
{
  return (p[2])/(sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]));
}

double phi(double p[ ])
{
  double phi;
  phi=atan2(p[1],p[0]);
  return phi;
}

void rotz(double p[],double theta)
{
  double px=cos(theta)*p[0]-sin(theta)*p[1];
  double py=sin(theta)*p[0]+cos(theta)*p[1];
  p[0] = px; p[1] = py;
  // p[2]=p[2];

}

void roty(double p[],double theta)
{
  double px= cos(theta)*p[0]+sin(theta)*p[2];
  double pz=-sin(theta)*p[0]+cos(theta)*p[2];
  p[0] = px; p[2] = pz;
  // p[2]=p[2];
  // return p;
}

void boost(double p[], double Gamma)
{
        double BG = Gamma > 0 ? sqrt(Gamma*Gamma-1.0) : -sqrt(Gamma*Gamma-1.0);
        double G = Gamma > 0 ? Gamma : -Gamma;
        double pz =  G*p[2] + BG*p[3];
  			double E  = BG*p[2] +  G*p[3];
        p[2] = pz;
        p[3] = E;
}

#define MASSPI 0.13957
#define MASSPROT 0.938
#define MASSPISQ (MASSPI*MASSPI)
#define MASSPROTSQ (MASSPROT*MASSPROT)

#define LAMBDA(x,y,z) ((x)*(x)+(y)*(y)+(z)*(z)-2*(x)*(y)-2*(y)*(z)-2*(z)*(x))

double tmin(double mXsq, double stot) {
  return MASSPISQ+mXsq-(stot+MASSPISQ-MASSPROTSQ)*(stot+mXsq-MASSPROTSQ)/(2*stot)+
    sqrt(LAMBDA(stot,MASSPISQ,MASSPROTSQ)*LAMBDA(stot,mXsq,MASSPROTSQ))/(2*stot);
}

double tprime(double pb[], double pX[]) {
  double mXsq = invmasssq(pX);
  double stot = MASSPISQ+MASSPROTSQ + 2*MASSPROTSQ*pb[3];
  double tmin_here = tmin(mXsq, stot);

  double diff[4]; for (uint i=0; i<4; i++) diff[i] = pb[i]-pX[i];
  double t = invmasssq(diff);  // negative

  return tmin_here-t;  // should be positive
}

void add(double p1[], double p2[], double sum[]) {
  for (int i = 0; i < 4; i++)
  sum[i] = p1[i] + p2[i];
}

void show(double p[]) {
  std::cout << "(px,py,pz,E) = ("
            << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << "), M = "
            << sqrt(invmasssq(p)) << "\n";
}

void apply_to_all(std::function<void(double*)> f, std::vector<double*> arr) {
  for (auto v : arr) f(v);
}


void print_ph_parameters(double p0[], double p1[], double p2[], double p3[],std::ofstream *BinFile) {

  double p123[4], p23[4];
  double pt[] = {0,0,0,MASSPROT};

  add(p2,p3,p23);
  add(p1,p23,p123);
  double s = invmasssq(p123);
  double sigma1 = invmasssq(p23);
  double th123 = acos(costheta(p123));
  double ph123 = phi(p123);
  apply_to_all([&](double *p)->void{rotz(p,-ph123); roty(p,-th123);}, {p1,p2,p3,p0,pt});

  // boost to CMS
  double gamma1=p123[3]/(sqrt(invmasssq(p123)));
  apply_to_all([&](double *p)->void{boost(p,-gamma1);},
               {p1,p2,p3,p0,pt});

  // rot to align beam along z
  double thb = acos(costheta(p0));
  double phb = phi(p0);
  apply_to_all([&](double *p)->void{rotz(p,-phb); roty(p,-thb); rotz(p,M_PI);},
               {p1,p2,p3,p0,pt});

  // measure Omega1
  add(p2,p3,p23);
  double th1=acos(costheta(p23));
  double theta1=costheta(p23);
  double ph1=phi(p23);
  // rot to align (23) to z axis
  apply_to_all([&](double *p)->void{rotz(p,-ph1); roty(p,-th1);},
               {p1,p2,p3,p0});
  // boost to 2,3 frame
  add(p2,p3,p23);
  double gamma23=p23[3]/(sqrt(invmasssq(p23)));
  apply_to_all([&](double *p)->void{boost(p,-gamma23);}, {p1,p2,p3,p0});
  // measure Omega23
  double th23=costheta(p2);
  double ph23=phi(p2);


  //std::cout << s << " " << sigma1 << " " << theta1 << " " << ph1 << " " << th23 << " " << ph23 << "\n";
  // std::cout << s << sigma1 << theta1 << ph1 << th23 << ph23 <<;
  BinFile->write((char *)&s, sizeof(double));
  BinFile->write((char *)&sigma1, sizeof(double));
  BinFile->write((char *)&theta1, sizeof(double));
  BinFile->write((char *)&ph1, sizeof(double));
  BinFile->write((char *)&th23, sizeof(double));
  BinFile->write((char *)&ph23, sizeof(double));
}

int ExtractMC(const char* file, bool full, double tprime_cut_L, double tprime_cut_H, std::ofstream * BinFile) {
  double p0[4],p1[4],p2[4],p3[4];
  TFile *fin = TFile::Open(file); if (!fin) return 1;
  TTree *tree = (TTree*)fin->Get("USR51MCout"); if (!tree) return 1;


  long nentries = tree->GetEntries();

  tree->SetBranchAddress("px0_MCTruth", &p0[0]);
  tree->SetBranchAddress("py0_MCTruth", &p0[1]);
  tree->SetBranchAddress("pz0_MCTruth", &p0[2]);
  tree->SetBranchAddress("px1_MCTruth", &p1[0]);
  tree->SetBranchAddress("py1_MCTruth", &p1[1]);
  tree->SetBranchAddress("pz1_MCTruth", &p1[2]);
  tree->SetBranchAddress("px2_MCTruth", &p2[0]);
  tree->SetBranchAddress("py2_MCTruth", &p2[1]);
  tree->SetBranchAddress("pz2_MCTruth", &p2[2]);
  tree->SetBranchAddress("px3_MCTruth", &p3[0]);
  tree->SetBranchAddress("py3_MCTruth", &p3[1]);
  tree->SetBranchAddress("pz3_MCTruth", &p3[2]);

  // cuts
  int accepted_reco;

  tree->SetBranchAddress("accepted_reco", &accepted_reco);
  int l=0;
  for (int i=0;i<nentries;i++) {
    static bool first=true;
    static bool second =false;
    tree->GetEntry(i);
    if (!full && accepted_reco != 1){
      if((i+1)==nentries){
      first=false;
      if(second) continue;
      BinFile->write((char *)&l, sizeof(int));
      second=true;
      i=0;
    }continue;}

    // process with event

    p0[3] = energy(p0,MASSPISQ);
    p1[3] = energy(p1,MASSPISQ);
    p2[3] = energy(p2,MASSPISQ);
    p3[3] = energy(p3,MASSPISQ);

    double pX[3]; add(p1,p2,pX); add(p3,pX,pX);
    double tprime_here = tprime(p0, pX);
    if (tprime_here < tprime_cut_L || tprime_here > tprime_cut_H){
      if((i+1)==nentries){
      first=false;
      if(second) continue;
      BinFile->write((char *)&l, sizeof(int));
      second=true;
      i=0;
    }continue;}

    if(first){if((i+1)==nentries){
    first=false;
    BinFile->write((char *)&l, sizeof(int));
    i=0;
    }l++;}
    else print_ph_parameters(p0,p2,p1,p3,BinFile);
  }

  // for (int i=0;i<nentries;i++) {
  //   tree->GetEntry(i);
  //   if (!full && accepted_reco != 1) continue;
  //
  //   // process with event
  //
  //   p0[3] = energy(p0,MASSPISQ);
  //   p1[3] = energy(p1,MASSPISQ);
  //   p2[3] = energy(p2,MASSPISQ);
  //   p3[3] = energy(p3,MASSPISQ);
  //
  //   double pX[3]; add(p1,p2,pX); add(p3,pX,pX);
  //   double tprime_here = tprime(p0, pX);
  //   if (tprime_here < tprime_cut_L || tprime_here > tprime_cut_H) continue;
  //
  //   print_ph_parameters(p0,p2,p1,p3,BinFile);
  //
  // }
  return 0;

}


int ExtractRD(const char* file, double tprime_cut_L, double tprime_cut_H,std::ofstream *BinFile) {
  double p0[4],p1[4],p2[4],p3[4];
  TFile *fin = TFile::Open(file); if (!fin) return 1;
  TTree *tree = (TTree*)fin->Get("USR52mb"); if (!tree) return 1;

  int k=0;
  long nentries = tree->GetEntries();

  tree->SetBranchAddress("px0", &p0[0]);
  tree->SetBranchAddress("py0", &p0[1]);
  tree->SetBranchAddress("pz0", &p0[2]);
  tree->SetBranchAddress("px1", &p1[0]);
  tree->SetBranchAddress("py1", &p1[1]);
  tree->SetBranchAddress("pz1", &p1[2]);
  tree->SetBranchAddress("px2", &p2[0]);
  tree->SetBranchAddress("py2", &p2[1]);
  tree->SetBranchAddress("pz2", &p2[2]);
  tree->SetBranchAddress("px3", &p3[0]);
  tree->SetBranchAddress("py3", &p3[1]);
  tree->SetBranchAddress("pz3", &p3[2]);

  // cuts
  bool IsTriggered, IsInTarget, IsExclusive, IsInT,
    CentralProdVeto, IsInBeamTime, RICH_Veto, CEDAR_Veto, IsInDeltaPhi,
    CorrectNbrRPDTracks, IsPlanar, IsPlanar_extended;

  tree->SetBranchAddress("IsTriggered",         &IsTriggered);
  tree->SetBranchAddress("IsInTarget",          &IsInTarget);
  tree->SetBranchAddress("IsExclusive",         &IsExclusive);
  tree->SetBranchAddress("IsInT" ,              &IsInT);
  tree->SetBranchAddress("CentralProdVeto",     &CentralProdVeto);
  tree->SetBranchAddress("IsInBeamTime",        &IsInBeamTime);
  tree->SetBranchAddress("RICH_Veto",           &RICH_Veto);
  tree->SetBranchAddress("CEDAR_Veto",          &CEDAR_Veto);
  tree->SetBranchAddress("IsInDeltaPhi",        &IsInDeltaPhi);
  tree->SetBranchAddress("CorrectNbrRPDTracks", &CorrectNbrRPDTracks);
  tree->SetBranchAddress("IsPlanar",            &IsPlanar);
  tree->SetBranchAddress("IsPlanar_extended",   &IsPlanar_extended);

  for (int i=0;i<nentries;i++) {
    static bool first = true;
    static bool second = false;
    tree->GetEntry(i);
    if (
        !IsTriggered ||
        !IsInTarget  ||
        !IsExclusive ||
        !IsInT  ||
        CentralProdVeto ||
        !IsInBeamTime ||
        RICH_Veto ||
        CEDAR_Veto ||
        // IsInDeltaPhi != 0 ||
        !CorrectNbrRPDTracks ||
        // IsPlanar != 0 ||
        !IsPlanar_extended
        ){
          if((i+1) == nentries){
          first=false;
          if(second) continue;
          BinFile->write((char *)&k, sizeof(int));
          second=true;
          i=0;
        }continue;}

    // process with event

    p0[3] = energy(p0,MASSPISQ);
    p1[3] = energy(p1,MASSPISQ);
    p2[3] = energy(p2,MASSPISQ);
    p3[3] = energy(p3,MASSPISQ);

    double pX[3]; add(p1,p2,pX); add(p3,pX,pX);
    double tprime_here = tprime(p0, pX);

    if (tprime_here < tprime_cut_L || tprime_here > tprime_cut_H) {
      if((i+1)==nentries) {
      first =false;
      if(second) continue;
      BinFile->write((char *)&k, sizeof(int));
      second=true;
      i=0;
    }continue;}

    if(first){
      if((i+1) == nentries){
        first=false;
        BinFile->write((char *)&k, sizeof(int));
        i=0;
    }k++;}
    else print_ph_parameters(p0,p2,p3,p1,BinFile);

  }

  return 0;
}

int Extract(const char* file, char flag,
            double tprime_cut_L, double tprime_cut_H, std::ofstream *BinFile) {
  if (flag != 'D' && flag != 'M' && flag != 'F') {
    std::cout << "There are two options for the second argument"
              << "\n\t'D' for real data"
              << "\n\t'M' for the accepted monteCarlo data"
              << "\n\t'M' for the full phase space monteCarlo data"
              << "\nYou gave --" << flag << "--\n";
    return 0;
  }
  if (flag == 'D') ExtractRD(file,        tprime_cut_L, tprime_cut_H,BinFile);
  if (flag == 'M') ExtractMC(file, false, tprime_cut_L, tprime_cut_H,BinFile);
  if (flag == 'F') ExtractMC(file, true,  tprime_cut_L, tprime_cut_H,BinFile);
  return 0;
}

int main(int argv, char **argc) {
  std::ofstream BinFile;
  BinFile.open(argc[5],std::ios::out | std::ios::binary);
  if (argv < 3) {
    std::cerr << "Arguments are needed!\n"
             << "Example:\n\t ./Extract datafile.root 'D'\n";
    return 1.0;
  }

  double tprime_cut_L = (argv < 5) ? 0.0 : atof(argc[3]);
  double tprime_cut_H = (argv < 5) ? 1.0 : atof(argc[4]);

  const char* file = argc[1];
  char flag = argc[2][0];

  Extract(file, flag, tprime_cut_L, tprime_cut_H,&BinFile);

  BinFile.close();

  return 0.0;
}
