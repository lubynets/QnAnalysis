#include "ShapeFitter.hpp"

ShapeFitter::ShapeFitter(TH1F* histo)
{
  histo_all_ = histo;
}

void ShapeFitter::Fit()
{
  const float mu = 1.11572;
  const float sigma = 0.00145786;
  
  const float left_external = mu - 13*sigma;
  const float left_internal = mu - 7*sigma;
  const float right_internal = mu + 7*sigma;
  const float right_external = mu + 13*sigma;

  TH1F* histo_without_peak = ExcludeInterval(histo_all_, left_internal, right_internal);
  TF1* func_bckgr = FitBckgr(histo_without_peak, left_external, right_external);
  func_bckgr_ = func_bckgr;
  histo_sgnl_ = SubtractBckgr(histo_all_, func_bckgr);
}

TH1F* ShapeFitter::ExcludeInterval(TH1F* histo, float left, float right)
{
  TH1F* histo_out = (TH1F*)histo->Clone();
  for(int iBin=histo_out->FindBin(left); iBin<=histo_out->FindBin(right); iBin++)
    histo_out -> SetBinContent(iBin, 0);
  
  return histo_out;
}

TF1* ShapeFitter::FitBckgr(TH1F* histo, float left, float right)
{
  TF1* bckgr_fit = new TF1 ("bckgr_fit", "pol2", left, right);
  histo -> Fit(bckgr_fit, "R0");
  
  return bckgr_fit;
}

TH1F* ShapeFitter::SubtractBckgr(TH1F* histo, TF1* func)
{
  TH1F* histo_out = (TH1F*)histo->Clone();
  histo_out -> Sumw2();
  histo_out -> Add(func, -1);
  
  return histo_out;
}