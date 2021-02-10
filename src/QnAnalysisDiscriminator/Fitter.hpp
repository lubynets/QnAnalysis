#ifndef Fitter_H
#define Fitter_H

#include"TH1.h"
#include"TFile.h"
#include"TString.h"
#include"TGraph.h"
#include"TF1.h"

class MyFunctor
{
public:  
  
  MyFunctor(TH1F* h1, TH1F* h2)
  {
    alpha_ = h1;
    beta_ = h2;
  }
  
  virtual ~MyFunctor() = default;
  
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
    return (GetTH1Value(alpha_, x[0])*signal_fit(x, par) + GetTH1Value(beta_, x[0])*bckgr_fit(x, &par[1])) / (GetTH1Value(alpha_, x[0]) + GetTH1Value(beta_, x[0]));
  }
  
private:
  
  const float mu = 1.11572;   //TODO remove this ugly hardcode ASAP
  
  double GetTH1Value(TH1F* h, double arg)
  {
    return h->GetBinContent(h->FindBin(arg));
  }
  
  TH1F* alpha_{nullptr};
  TH1F* beta_{nullptr};
};

class Fitter
{
public:
  
  Fitter() = default;
  virtual ~Fitter() = default;  
  
  void SetSignalShape(TH1F* histo) { histo_sgnl_ = histo; };
  void SetBackgroundShape(TH1F* histo) { histo_bckgr_ = histo; };
  void SetGraphToFit(TGraph* graph) { graph_v_ = graph; };
  
  void Fit();
  
  
private:
  
  TH1F* histo_sgnl_{nullptr};
  TH1F* histo_bckgr_{nullptr};
  TGraph* graph_v_{nullptr};
};

#endif//Fitter_H