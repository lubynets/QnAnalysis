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

// TF1* ShapeFitter::FitSgnl(TH1F* histo, float left, float right) const // try to fit with polynom8. Pre-definition of coefficients by solving SLE
// {  
//   TMatrixD A(9, 9);
//   
//   std::array<float, 9> x;
//   TVectorD b(9);
//   TDecompLU lu(9);
//   
//   x.at(0) = 0. - 14*ShapeFitter::sigma;
//   x.at(1) = 0. - 10*ShapeFitter::sigma;
//   x.at(2) = 0. -  6*ShapeFitter::sigma;
//   x.at(3) = 0. -  3*ShapeFitter::sigma;
//   x.at(4) = 0.                        ;
//   x.at(5) = 0. +  3*ShapeFitter::sigma;
//   x.at(6) = 0. +  6*ShapeFitter::sigma;
//   x.at(7) = 0. + 10*ShapeFitter::sigma;
//   x.at(8) = 0. + 14*ShapeFitter::sigma;
//   
//   b[0] = histo->Interpolate(ShapeFitter::mu - 14*ShapeFitter::sigma);
//   b[1] = histo->Interpolate(ShapeFitter::mu - 10*ShapeFitter::sigma);
//   b[2] = histo->Interpolate(ShapeFitter::mu -  6*ShapeFitter::sigma);
//   b[3] = histo->Interpolate(ShapeFitter::mu -  3*ShapeFitter::sigma);
//   b[4] = histo->Interpolate(ShapeFitter::mu                        );
//   b[5] = histo->Interpolate(ShapeFitter::mu +  3*ShapeFitter::sigma);
//   b[6] = histo->Interpolate(ShapeFitter::mu +  6*ShapeFitter::sigma);
//   b[7] = histo->Interpolate(ShapeFitter::mu + 10*ShapeFitter::sigma);
//   b[8] = histo->Interpolate(ShapeFitter::mu + 14*ShapeFitter::sigma);
//   
//   for(int iPosition=0; iPosition<9; iPosition++)
//     for(int iPower=0; iPower<9; iPower++)
//       A[iPosition][iPower] = std::pow(x.at(iPosition), iPower);
//   
//   A.Invert();
//   
//   TVectorD p = A*b;
//   
//   p.Print();
//   
//   TF1* sgnl_fit = new TF1("sgnl_fit", "[0]+[1]*(x-[9])+[2]*(x-[9])*(x-[9])+[3]*(x-[9])**3+[4]*(x-[9])**4+[5]*(x-[9])**5+[6]*(x-[9])**6+[7]*(x-[9])**7+[8]*(x-[9])**8", left, right);
//   
//   sgnl_fit -> FixParameter(9, ShapeFitter::mu);
//   
//   for(int i=0; i<9; i++)
//     sgnl_fit -> SetParameter(i, p[i]);
//   
//   return sgnl_fit;
// }


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


// //********** two expos and gauss ******************************************************
// TF1* ShapeFitter::FitSgnl(TH1F* histo, float left, float right) const
// {
//   const int Npar = 6;
//   
//   MyFunctorShape myfuncsh;
//   TF1* sgnl_fit = new TF1("sgnl_fit", myfuncsh, left, right, Npar);
//   sgnl_fit -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::Gaus(0, 0, ShapeFitter::sigma));
//   sgnl_fit -> FixParameter(1, ShapeFitter::mu);
//   sgnl_fit -> SetParameter(2, 0);
//   sgnl_fit -> SetParameter(3, ShapeFitter::sigma);
//   sgnl_fit -> SetParameter(4, 1000.);
//   sgnl_fit -> SetParameter(5, -1000.);
//   
//   histo -> Fit(sgnl_fit, "R0");
// 
//   return sgnl_fit;
// }
// 
// double MyFunctorShape::operator()(double* x, double* par)
// {
//   const double factor     = par[0];
//   const double shift      = par[1];
//   const double mu         = par[2]; // expected to be 0
//   const double sigma      = par[3];
//   const double k_left     = par[4];
//   const double k_right    = par[5];
//   
//   const double xx = x[0] - shift;
//     
//   const double C_left  = TMath::Exp(-x0_*x0_/2/sigma/sigma + k_left*x0_);
//   const double C_right = TMath::Exp(-x0_*x0_/2/sigma/sigma - k_right*x0_);
//   
//   if(xx < -x0_)
//     return factor*C_left*TMath::Exp(k_left * xx);
//   else if (xx > x0_)
//     return factor*C_right*TMath::Exp(k_right * xx);
//   else
//     return factor*TMath::Gaus(xx, mu, sigma);  
// }
// //*************************************************************************************

//********** two expos and lorentz ****************************************************
TF1* ShapeFitter::FitSgnl(TH1F* histo, float left, float right) const
{
  const int Npar = 6;
  
  MyFunctorShape myfuncsh;
  TF1* sgnl_fit = new TF1("sgnl_fit", myfuncsh, left, right, Npar);
  sgnl_fit -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::CauchyDist(0, 0, ShapeFitter::sigma));
  sgnl_fit -> FixParameter(1, ShapeFitter::mu);
  sgnl_fit -> SetParameter(2, 0);
  sgnl_fit -> SetParameter(3, ShapeFitter::sigma);
  sgnl_fit -> SetParameter(4, 1000.);
  sgnl_fit -> SetParameter(5, -1000.);
  
  histo -> Fit(sgnl_fit, "R0");

  return sgnl_fit;
}

double MyFunctorShape::operator()(double* x, double* par)
{
  const double factor     = par[0];
  const double shift      = par[1];
  const double mu         = par[2]; // expected to be 0
  const double sigma      = par[3];
  const double k_left     = par[4];
  const double k_right    = par[5];
  
  const double xx = x[0] - shift;
    
  const double C_left  = TMath::CauchyDist(-x0_, 0, sigma) / TMath::Exp(-k_left*x0_);
  const double C_right = TMath::CauchyDist( x0_, 0, sigma) / TMath::Exp( k_right*x0_);
  
  if(xx < -x0_)
    return factor*C_left*TMath::Exp(k_left * xx);
  else if (xx > x0_)
    return factor*C_right*TMath::Exp(k_right * xx);
  else
    return factor*TMath::CauchyDist(xx, mu, sigma);  
}
//*************************************************************************************