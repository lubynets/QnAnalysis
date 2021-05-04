#include "ShapeFitter.hpp"

ShapeFitter::ShapeFitter(TH1F* histo)
{
  histo_all_ = histo;
}

void ShapeFitter::Fit()
{
  const float left_external = ShapeFitter::mu - 15 * ShapeFitter::sigma;      //TODO make it settable
  const float left_internal = ShapeFitter::mu - 8 * ShapeFitter::sigma;
  const float right_internal = ShapeFitter::mu + 8 * ShapeFitter::sigma;
  const float right_external = ShapeFitter::mu + 15 * ShapeFitter::sigma;

  TH1F* histo_without_peak = ExcludeInterval(histo_all_, left_internal, right_internal);
  func_bckgr_ = FitBckgr(histo_without_peak, left_external, right_external);
  histo_sgnl_ = SubtractBckgr(histo_all_, func_bckgr_, left_external, right_external);
  chi2_bckgr_fit_ = func_bckgr_->GetChisquare() / func_bckgr_->GetNDF();
  
  func_sgnl_ = FitSgnl(histo_sgnl_, left_external, right_external);
  chi2_sgnl_fit_ = func_sgnl_->GetChisquare() / func_sgnl_->GetNDF();
}

TH1F* ShapeFitter::ExcludeInterval(TH1F* histo, float left, float right) const
{
  TH1F* histo_out = (TH1F*)histo->Clone();
  for(int iBin=histo_out->FindBin(left); iBin<=histo_out->FindBin(right); iBin++)
    histo_out -> SetBinContent(iBin, 0);
  
  return histo_out;
}

TF1* ShapeFitter::FitBckgr(TH1F* histo, float left, float right) const
{
  TF1* bckgr_fit = new TF1 ("bckgr_fit", "pol4", left, right);                            // TODO make it settable
  histo -> Fit(bckgr_fit, "R0");
    
  return bckgr_fit;
}

TF1* ShapeFitter::FitSgnl(TH1F* histo, float left, float right) const
{
  TF1* sgnl_fit = new TF1("sgnl_fit", "[0]*TMath::Gaus(x-[1], [2], [3])", left, right);
  sgnl_fit -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::Gaus(0, 0, ShapeFitter::sigma));
  sgnl_fit -> FixParameter(1, ShapeFitter::mu);
  sgnl_fit -> SetParameter(2, 0.);
  sgnl_fit -> SetParameter(3, ShapeFitter::sigma);  
  
//   TF1* sgnl_fit = new TF1("sgnl_fit", "[0]*TMath::CauchyDist(x-[1], [2], [3])", left, right);
//   sgnl_fit -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::CauchyDist(0, 0, ShapeFitter::sigma));
//   sgnl_fit -> FixParameter(1, ShapeFitter::mu);
//   sgnl_fit -> SetParameter(2, 0.);
//   sgnl_fit -> SetParameter(3, ShapeFitter::sigma);
//   
//   TF1* sgnl_fit = new TF1("sgnl_fit", "[0]*TMath::Voigt(x-[1]-[2], [3], [4])", left, right);
//   sgnl_fit -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::Voigt(0, ShapeFitter::sigma, ShapeFitter::sigma));
//   sgnl_fit -> FixParameter(1, ShapeFitter::mu);
//   sgnl_fit -> SetParameter(1, 0);
//   sgnl_fit -> SetParameter(2, ShapeFitter::sigma);
//   sgnl_fit -> SetParameter(3, ShapeFitter::sigma);
  
  histo -> Fit(sgnl_fit, "R0");
  
  return sgnl_fit;
}

TH1F* ShapeFitter::SubtractBckgr(TH1F* histo, TF1* func, float left, float right) const
{
  TH1F* histo_out = (TH1F*)histo->Clone();
  histo_out -> Sumw2();
  histo_out -> Add(func, -1);
  
  for(int iBin=1; iBin<=histo_out->GetNbinsX(); iBin++)
    if(histo_out->GetBinCenter(iBin)<left || histo_out->GetBinCenter(iBin)>right)
    {
      histo_out->SetBinContent(iBin, 0);
      histo_out->SetBinError(iBin, 0);
    }
      
  return histo_out;
}
