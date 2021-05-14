#include "ShapeFitter.hpp"

#include "TMatrixD.h"
#include "TFitResult.h"

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
  graph_bckgr_ = FuncWithErrors(func_bckgr_);
  chi2_bckgr_fit_ = func_bckgr_.first->GetChisquare() / func_bckgr_.first->GetNDF();
  
  histo_sgnl_ = SubtractBckgr(histo_all_, func_bckgr_, left_external, right_external);  
  func_sgnl_ = FitSgnl(histo_sgnl_, left_external, right_external);
  graph_sgnl_ = FuncWithErrors(func_sgnl_);
  chi2_sgnl_fit_ = func_sgnl_.first->GetChisquare() / func_sgnl_.first->GetNDF();
}

TH1F* ShapeFitter::ExcludeInterval(TH1F* histo, float left, float right) const
{
  TH1F* histo_out = (TH1F*)histo->Clone();
  for(int iBin=histo_out->FindBin(left); iBin<=histo_out->FindBin(right); iBin++)
    histo_out -> SetBinContent(iBin, 0);
  
  return histo_out;
}

TH1F* ShapeFitter::SubtractBckgr(TH1F* histo, std::pair<TF1*, TMatrixDSym*> f_and_cov, float left, float right) const
{
  TH1F* histo_out = (TH1F*)histo->Clone();
  histo_out -> Sumw2();
  histo_out -> Add(f_and_cov.first, -1);
  
  for(int iBin=1; iBin<=histo_out->GetNbinsX(); iBin++)
  {
    if(histo_out->GetBinCenter(iBin)<left || histo_out->GetBinCenter(iBin)>right)
    {
      histo_out->SetBinContent(iBin, 0);
      histo_out->SetBinError(iBin, 0);
    }
    else
    {
      const float eh = histo_out->GetBinError(iBin);
      double binCenter = histo_out->GetBinCenter(iBin);
      const float ef = EvalError(&binCenter, f_and_cov);
      const float ee = std::sqrt(eh*eh + ef*ef);
      histo_out -> SetBinError(iBin, ee);
    }
  }
      
  return histo_out;
}

float ShapeFitter::EvalError(double* x, std::pair<TF1*, TMatrixDSym*> f_and_cov) const            // add check if npar of func is equal to dim cov
{
  const int Npar = f_and_cov.first->GetNpar();
  TMatrixD dfdp(Npar, 1);
  for(int i=0; i<Npar; i++)
    dfdp[i][0] = f_and_cov.first->GradientPar(i, x);
  
  TMatrixD dfdp_T = dfdp;
  dfdp_T.T();
  
  return std::sqrt((dfdp_T*(*f_and_cov.second)*dfdp)[0][0]);
}

TGraphErrors* ShapeFitter::FuncWithErrors(std::pair<TF1*, TMatrixDSym*> f_and_cov) const
{
  TGraphErrors* graph = new TGraphErrors();
  const int Nsteps = 1000;
  const float left = f_and_cov.first->GetXmin();
  const float right = f_and_cov.first->GetXmax();
  const float step = (right-left)/Nsteps;
  
  double x = left;
  int i = 0;
  while(x<=right)
  {
    const float y = f_and_cov.first->Eval(x);
    const float ey = EvalError(&x, f_and_cov);
    graph->SetPoint(i, x, y);
    graph->SetPointError(i, 0, ey);
    x += step;
    i++;
  }  
  
  return graph;
}

std::pair<TF1*, TMatrixDSym*>  ShapeFitter::FitBckgr(TH1F* histo, float left, float right) const
{
  TF1* bckgr_fit = new TF1("bckgr_fit", "pol2", left, right);                            // TODO make it settable
  TMatrixDSym* cov = new TMatrixDSym(bckgr_fit->GetNpar());
  
  TFitResultPtr frptr = histo -> Fit(bckgr_fit, "RS0");
  *cov = frptr -> GetCovarianceMatrix();
  
  std::pair<TF1*, TMatrixDSym*> f_and_cov(bckgr_fit, cov);
  
  return f_and_cov;
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
std::pair<TF1*, TMatrixDSym*> ShapeFitter::FitSgnl(TH1F* histo, float left, float right) const
{
  const int Npar = 6;
  
  MyFunctorShape myfuncsh;
  TF1* sgnl_fit = new TF1("sgnl_fit", myfuncsh, left, right, Npar);
  sgnl_fit -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::Gaus(0, 0, ShapeFitter::sigma));
  sgnl_fit -> FixParameter(1, ShapeFitter::mu);
  sgnl_fit -> SetParameter(2, 0);
  sgnl_fit -> SetParameter(3, ShapeFitter::sigma);
  sgnl_fit -> SetParameter(4, 1000.);
  sgnl_fit -> SetParameter(5, -1000.);
  
  TMatrixDSym* cov = new TMatrixDSym(sgnl_fit->GetNpar());
  
  TFitResultPtr frptr = histo -> Fit(sgnl_fit, "RS0");
  *cov = frptr -> GetCovarianceMatrix();
  
  std::pair<TF1*, TMatrixDSym*> f_and_cov(sgnl_fit, cov);
  
  return f_and_cov;
}

double MyFunctorShape::operator()(double* x, double* par)
{
  const double factor     = par[0];
  const double shift      = par[1]; // to be fixed at real peak position
  const double mu         = par[2]; // expected to be 0
  const double sigma      = par[3];
  const double k_left     = par[4];
  const double k_right    = par[5];
  
  const double xx = x[0] - shift;
    
  const double C_left  = TMath::Gaus(-x0_, mu, sigma) / TMath::Exp(-k_left*x0_);
  const double C_right = TMath::Gaus( x0_, mu, sigma) / TMath::Exp( k_right*x0_);
  
  if(xx < -x0_)
    return factor*C_left*TMath::Exp(k_left * xx);
  else if (xx > x0_)
    return factor*C_right*TMath::Exp(k_right * xx);
  else
    return factor*TMath::Gaus(xx, mu, sigma);  
}
//*************************************************************************************

// //********** two expos and lorentz ****************************************************
// TF1* ShapeFitter::FitSgnl(TH1F* histo, float left, float right) const
// {
//   const int Npar = 6;
//   
//   MyFunctorShape myfuncsh;
//   TF1* sgnl_fit = new TF1("sgnl_fit", myfuncsh, left, right, Npar);
//   sgnl_fit -> SetParameter(0, histo->Interpolate(ShapeFitter::mu) / TMath::CauchyDist(0, 0, ShapeFitter::sigma));
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
//   const double shift      = par[1]; // to be fixed at real peak position
//   const double mu         = par[2]; // expected to be 0
//   const double sigma      = par[3];
//   const double k_left     = par[4];
//   const double k_right    = par[5];
//   
//   const double xx = x[0] - shift;
//     
//   const double C_left  = TMath::CauchyDist(-x0_, mu, sigma) / TMath::Exp(-k_left*x0_);
//   const double C_right = TMath::CauchyDist( x0_, mu, sigma) / TMath::Exp( k_right*x0_);
//   
//   if(xx < -x0_)
//     return factor*C_left*TMath::Exp(k_left * xx);
//   else if (xx > x0_)
//     return factor*C_right*TMath::Exp(k_right * xx);
//   else
//     return factor*TMath::CauchyDist(xx, mu, sigma);  
// }
// //*************************************************************************************