#include "Fitter.hpp"

#include "TFile.h"
#include "TCanvas.h"

void Fitter::Fit()
{
  const float graphleft = graph_v_ -> GetPointX(0);                     //TODO read this info from Qn::Axis
  const float graphright = graph_v_ -> GetPointX(graph_v_->GetN()-1);

  MyFunctor funct(shape_);
  TF1* f = new TF1("f", funct, graphleft-0.01, graphright+0.01, 3);
  
//   f->SetParLimits(0, -2, 2);
//   f->SetParLimits(1, -3, 3);
//   f->SetParLimits(2, -20, 20);
//   f->SetParLimits(3, -50, 50);
    
  graph_v_ -> Fit("f");
  
  for(int i=0; i<3; i++)
  {
    fit_params_.push_back(f->GetParameter(i));
    fit_params_errors_.push_back(f->GetParError(i));
  }
}