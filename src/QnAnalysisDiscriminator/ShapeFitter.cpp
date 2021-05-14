#include "ShapeFitter.hpp"

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompLU.h"

#include <array>

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
  
  func_sgnl_ = FitSgnl(histo_sgnl_, left_external, right_external);        //TODO uncomment after debug
  chi2_sgnl_fit_ = func_sgnl_->GetChisquare() / func_sgnl_->GetNDF();
}

TH1F* ShapeFitter::ExcludeInterval(TH1F* histo, float left, float right) const
{
  TH1F* histo_out = (TH1F*)histo->Clone();
  for(int iBin=histo_out->FindBin(left); iBin<=histo_out->FindBin(right); iBin++)
    histo_out -> SetBinContent(iBin, 0);
  
  return histo_out;
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

TF1* ShapeFitter::FitBckgr(TH1F* histo, float left, float right) const
{
  TF1* bckgr_fit = new TF1 ("bckgr_fit", "pol3", left, right);                            // TODO make it settable
  histo -> Fit(bckgr_fit, "R0");
    
  return bckgr_fit;
}

// TF1* ShapeFitter::FitSgnl(TH1F* histo, float left, float right) const
// {
//   TF1* sgnl_fit = new TF1("sgnl_fit", "[0]*TMath::Gaus(x-[1], [2], [3])", left, right);
//   sgnl_fit -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::Gaus(0, 0, ShapeFitter::sigma));
//   sgnl_fit -> FixParameter(1, ShapeFitter::mu);
//   sgnl_fit -> SetParameter(2, 0.);
//   sgnl_fit -> SetParameter(3, ShapeFitter::sigma);  
//   
//   TF1* sgnl_fit = new TF1("sgnl_fit", "[0]*TMath::CauchyDist(x-[1], [2], [3])", left, right);
//   sgnl_fit -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::CauchyDist(0, 0, ShapeFitter::sigma));
//   sgnl_fit -> FixParameter(1, ShapeFitter::mu);
//   sgnl_fit -> SetParameter(2, 0.);
//   sgnl_fit -> SetParameter(3, ShapeFitter::sigma);
//    
//   TF1* sgnl_fit = new TF1("sgnl_fit", "[0]*TMath::Voigt(x-[1]-[2], [3], [4])", left, right);
//   sgnl_fit -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::Voigt(0, ShapeFitter::sigma/2, ShapeFitter::sigma/2));
//   sgnl_fit -> FixParameter(1, ShapeFitter::mu);
//   sgnl_fit -> SetParameter(2, 0);
//   sgnl_fit -> SetParameter(3, ShapeFitter::sigma/2);
//   sgnl_fit -> SetParameter(4, ShapeFitter::sigma/2);
//   
//   histo -> Fit(sgnl_fit, "R0");
//   
//   return sgnl_fit;
// }  

//********** two expos and gauss ******************************************************
TF1* ShapeFitter::FitSgnl(TH1F* histo, float left, float right) const
{
  const int Npar = 8;
  
  MyFunctorShape myfuncsh;
  TF1* sgnl_fit = new TF1("sgnl_fit", myfuncsh, left, right, Npar);
  sgnl_fit -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::Gaus(0, 0, ShapeFitter::sigma));
  sgnl_fit -> FixParameter(1, ShapeFitter::mu);
  sgnl_fit -> SetParameter(2, 0);
  sgnl_fit -> SetParameter(3, ShapeFitter::sigma);
  const float x_shift = 1.3e-3;     // point where Exp() is pre-defined
  sgnl_fit -> SetParameter(4, histo->Interpolate(ShapeFitter::mu - x_shift) / TMath::Exp(-1000*x_shift));
  sgnl_fit -> SetParameter(5, histo->Interpolate(ShapeFitter::mu + x_shift) / TMath::Exp(-1000*x_shift));
  sgnl_fit -> SetParameter(6, 1000.);
  sgnl_fit -> SetParameter(7, -1000.);
  
  histo -> Fit(sgnl_fit, "R0");

  return sgnl_fit;
}

double MyFunctorShape::operator()(double* x, double* par)
{
  const double factor_peak  = par[0];
  const double shift        = par[1]; // to be fixed at real peak position
  const double mu           = par[2]; // expected to be 0
  const double sigma        = par[3];
  const double factor_left  = par[4];
  const double factor_right = par[5];
  const double k_left       = par[6];
  const double k_right      = par[7];
  
  const double xx = x[0] - shift;
  
  auto alpha = [xx]           // fraction of Peak funkcion in transition region
  {
    double ksi = std::abs(xx);
    if(ksi<x0_internal_)
      return 1.;
    else if(ksi>x0_external_)
      return 0.;
    else
      return (x0_external_ - ksi) / (x0_external_ - x0_internal_);
  };
  
  if(xx < -x0_external_)
    return factor_left*TMath::Exp(k_left * xx);
  else if(xx > x0_external_)
    return factor_right*TMath::Exp(k_right * xx);
  else if(std::abs(xx) < x0_internal_)
    return factor_peak*TMath::Gaus(xx, mu, sigma);
  else if(xx>x0_internal_ && xx<x0_external_)
    return alpha()*factor_peak*TMath::Gaus(xx, mu, sigma) + (1 - alpha())*factor_right*TMath::Exp(k_right * xx);
  else if(xx<-x0_internal_ && xx>-x0_external_)
    return alpha()*factor_peak*TMath::Gaus(xx, mu, sigma) + (1 - alpha())*factor_left*TMath::Exp(k_left * xx);
  else
    return 0.;
}
//*************************************************************************************

// //********** two expos and lorentz ****************************************************
// TF1* ShapeFitter::FitSgnl(TH1F* histo, float left, float right) const
// {
//   const int Npar = 8;
//   
//   MyFunctorShape myfuncsh;
//   TF1* sgnl_fit = new TF1("sgnl_fit", myfuncsh, left, right, Npar);
//   sgnl_fit -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::CauchyDist(0, 0, ShapeFitter::sigma));
//   sgnl_fit -> FixParameter(1, ShapeFitter::mu);
//   sgnl_fit -> SetParameter(2, 0);
//   sgnl_fit -> SetParameter(3, ShapeFitter::sigma);
//   const float x_shift = 1.3e-3;     // point where Exp() is pre-defined
//   sgnl_fit -> SetParameter(4, histo->Interpolate(ShapeFitter::mu - x_shift) / TMath::Exp(-1000*x_shift));
//   sgnl_fit -> SetParameter(5, histo->Interpolate(ShapeFitter::mu + x_shift) / TMath::Exp(-1000*x_shift));
//   sgnl_fit -> SetParameter(6, 1000.);
//   sgnl_fit -> SetParameter(7, -1000.);
//   
//   histo -> Fit(sgnl_fit, "R0");
// 
//   return sgnl_fit;
// }
// 
// double MyFunctorShape::operator()(double* x, double* par)
// {
//   const double factor_peak  = par[0];
//   const double shift        = par[1]; // to be fixed at real peak position
//   const double mu           = par[2]; // expected to be 0
//   const double sigma        = par[3];
//   const double factor_left  = par[4];
//   const double factor_right = par[5];
//   const double k_left       = par[6];
//   const double k_right      = par[7];
//   
//   const double xx = x[0] - shift;
//   
//   auto alpha = [xx]           // fraction of Peak funkcion in transition region
//   {
//     double ksi = std::abs(xx);
//     if(ksi<x0_internal_)
//       return 1.;
//     else if(ksi>x0_external_)
//       return 0.;
//     else
//       return (x0_external_ - ksi) / (x0_external_ - x0_internal_);
//   };
//   
//   if(xx < -x0_external_)
//     return factor_left*TMath::Exp(k_left * xx);
//   else if(xx > x0_external_)
//     return factor_right*TMath::Exp(k_right * xx);
//   else if(std::abs(xx) < x0_internal_)
//     return factor_peak*TMath::CauchyDist(xx, mu, sigma);
//   else if(xx>x0_internal_ && xx<x0_external_)
//     return alpha()*factor_peak*TMath::CauchyDist(xx, mu, sigma) + (1 - alpha())*factor_right*TMath::Exp(k_right * xx);
//   else if(xx<-x0_internal_ && xx>-x0_external_)
//     return alpha()*factor_peak*TMath::CauchyDist(xx, mu, sigma) + (1 - alpha())*factor_left*TMath::Exp(k_left * xx);
//   else
//     return 0.;
// }
// //*************************************************************************************