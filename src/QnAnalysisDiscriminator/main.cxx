#include <iostream>
#include"TH1.h"
#include"TFile.h"
#include"TString.h"
#include"TGraph.h"
#include"TF1.h"

class MyFunctor
{
public:
  const float mu = 1.11572;
  
  MyFunctor(TH1F* h1, TH1F* h2)
  {
    alpha_ = h1;
    beta_ = h2;
  }
  
  double signal_fit(double* x, double* par)
  {
    return par[0];
  }
  
  double bckgr_fit(double* x, double* par)
  {
    return par[0] + par[1]*(x[0]-mu) + par[2]*(x[0]-mu)*(x[0]-mu);
  }
    
  double operator()(double* x, double* par)
  {
    return (GetValue(alpha_, x[0])*signal_fit(x, par) + GetValue(beta_, x[0])*bckgr_fit(x, &par[1])) / (GetValue(alpha_, x[0]) + GetValue(beta_, x[0]));
  }
  
private:
  
  double GetValue(TH1F* h, double arg)
  {
    return h->GetBinContent(h->FindBin(arg));
  }
  
  TH1F* alpha_{nullptr};
  TH1F* beta_{nullptr};
};

int main(int argc, char** argv)
{
  TString sgnlfilename="/home/user/cbmdir/working/massfit/pfsqa.sgnl_12.root";
  TString bckgrfilename="/home/user/cbmdir/working/massfit/pfsqa.bckgr.root";
  TString v1filename="/home/user/cbmdir/working/massfit/v1graph.toymc.fmd.all.root";
  
  TFile* v1file = TFile::Open(v1filename, "read");
  TGraph* v1graph = (TGraph*)v1file -> Get("v1_invmass");
  
  TFile* sgnlfile = TFile::Open(sgnlfilename, "read");
  TH1F* histosgnl = (TH1F*)sgnlfile -> Get("LambdaCandidates/LambdaCandidates_mass");
  TFile* bckgrfile = TFile::Open(bckgrfilename, "read");
  TH1F* histobckgr = (TH1F*)bckgrfile -> Get("LambdaCandidates/LambdaCandidates_mass");
  
  const float graphleft = v1graph -> GetPointX(0);
  const float graphright = v1graph -> GetPointX(v1graph->GetN()-1);

  MyFunctor funct(histosgnl, histobckgr);
  TF1* f = new TF1("f", funct, graphleft-0.1, graphright+0.1, 4);
  
  v1graph -> Fit("f");
  
  v1graph -> Draw("AP");
  f -> Draw("same");
  
  return 0;
}