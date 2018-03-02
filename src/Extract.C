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
#define MASSPISQ (MASSPI*MASSPI)

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

int Extract(const char* file)
{
  double p0[4],p1[4],p2[4],p3[4];
  TFile *fin = TFile::Open(file); if (!fin) return 1;
  TTree *tree = (TTree*)fin->Get("USR52mb"); if (!tree) return 1;


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

  for (int i=0;i<nentries;i++)
  {
      tree->GetEntry(i);
      p0[3] = energy(p0,MASSPISQ);
      p1[3] = energy(p1,MASSPISQ);
      p2[3] = energy(p2,MASSPISQ);
      p3[3] = energy(p3,MASSPISQ);


      double p123[4], p23 [4];


      add(p2,p3,p23);
      add(p1,p23,p123);
      double s = invmasssq(p123);
      double sigma1 = invmasssq(p23);
      double th123 = acos(costheta(p123));
      double ph123 = phi(p123);
      // rotz(p1,-ph123); roty(p1,-th123);

      apply_to_all([&](double *p)->void{rotz(p,-ph123); roty(p,-th123);}, {p1,p2,p3,p0});

      // boost to CMS
      double gamma1=p123[3]/(sqrt(invmasssq(p123)));
      apply_to_all([&](double *p)->void{boost(p,-gamma1);},
                   {p1,p2,p3,p0});
      // rot to align beam along z
      double thb = acos(costheta(p0));
      double phb = phi(p0);
      apply_to_all([&](double *p)->void{rotz(p,-phb); roty(p,-thb);},
                   {p1,p2,p3,p0});



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


         cout << s << " " << sigma1 << " " << theta1 << " " << ph1 << " " << th23 << " " << ph23 << "\n";

   }

  return 0;
}
