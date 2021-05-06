#include "Fitter.hpp"

#include "TFile.h"
#include "TCanvas.h"

void Fitter::Fit()
{
  const int Npar = 3;
  
  const float graphleft = graph_v_ -> GetPointX(0);                     //TODO read this info from Qn::Axis
  const float graphright = graph_v_ -> GetPointX(graph_v_->GetN()-1);

  MyFunctor funct(shape_);
  TF1* f = new TF1("f", funct, graphleft-0.01, graphright+0.01, Npar);
     
  graph_v_ -> Fit("f");
  fit_chi2_ = f->GetChisquare();
  fit_ndf_ = f->GetNDF();
  
  TF1* f2 = new TF1("f2", "[0] + [1]*(x-1.11572)", graphleft-0.01, graphright+0.01);  // contribution from bckgr to flow
  
  for(int i=0; i<Npar; i++)
  {
    if(i!=0)
      f2->SetParameter(i-1, f->GetParameter(i));
    
    fit_params_.push_back(f->GetParameter(i));
    fit_params_errors_.push_back(f->GetParError(i));
  }
  
  f2 -> SetLineColor(kBlue);
  graph_v_ -> GetListOfFunctions() -> Add(f2);
}