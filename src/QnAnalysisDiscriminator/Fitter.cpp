#include "Fitter.hpp"

#include "TFile.h"
#include "TCanvas.h"

void Fitter::Fit()
{
  const float graphleft = graph_v_ -> GetPointX(0);                     //TODO read this info from Qn::Axis
  const float graphright = graph_v_ -> GetPointX(graph_v_->GetN()-1);

  MyFunctor funct(shape_);
  TF1* f = new TF1("f", funct, graphleft-0.01, graphright+0.01, 4);
  
  graph_v_ -> Fit("f");
  
  TFile* fileOut = TFile::Open("fit.root", "recreate");
  TCanvas* cc = new TCanvas("cc", "CC", 900, 1200);
  cc -> cd();
  
  graph_v_ -> Draw("AP");
  f -> SetLineColor(kBlue);
  f -> Draw("same");
  
  cc -> Write();
  fileOut -> Close();
}