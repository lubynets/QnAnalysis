#include "Fitter.hpp"

void Fitter::Fit()
{
  const float graphleft = graph_v_ -> GetPointX(0);                     //TODO read this info from Qn::Axis
  const float graphright = graph_v_ -> GetPointX(graph_v_->GetN()-1);

  MyFunctor funct(histo_sgnl_, histo_bckgr_);
  TF1* f = new TF1("f", funct, graphleft-0.1, graphright+0.1, 4);
  
  graph_v_ -> Fit("f");
  
  graph_v_ -> Draw("AP");
  f -> Draw("same");  
}