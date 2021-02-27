#include "ShapeContainer.hpp"

ClassImp(ShapeContainer);

void ShapeContainer::SetShape(TH1F* histosgnl, TH1F* histobckgr)
{
  histo_sgnl_ = histosgnl;
  histo_bckgr_ = histobckgr;
  is_histo_ = true;
}

void ShapeContainer::SetShape(TF1* funcsgnl, TF1* funcbckgr)
{
  func_sgnl_ =  funcsgnl;
  func_bckgr_ = funcbckgr;
  is_histo_ = false;
}

float ShapeContainer::GetSignal(float x)                       //TODO add warnings or even errors if histo or/and function is not initialized
{
  if(is_histo_ == true)
    return histo_sgnl_ -> Interpolate(x);
  else
    return func_sgnl_ -> Eval(x);
}

float ShapeContainer::GetBackground(float x)                   //TODO add warnings or even errors if histo or/and function is not initialized
{
  if(is_histo_ == true)
    return histo_bckgr_ -> Interpolate(x);
  else
    return func_bckgr_ -> Eval(x);
}

float ShapeContainer::GetSignalIntegral(float left, float right)
{
  if(is_histo_ == true)
    return ShapeContainer::HistoIntegral(histo_sgnl_, left, right);
  else
    return func_sgnl_ -> Integral(left, right);
}

float ShapeContainer::GetBackgroundIntegral(float left, float right)
{
  if(is_histo_ == true)
    return ShapeContainer::HistoIntegral(histo_bckgr_, left, right);
  else
    return func_bckgr_ -> Integral(left, right);
}

float ShapeContainer::HistoIntegral(TH1F* histo, float low, float up)
{
  if(up<low)
    return -ShapeContainer::HistoIntegral(histo, up, low);
  
  if(low<histo->GetXaxis()->GetBinLowEdge(0))
    low = histo->GetXaxis()->GetBinLowEdge(0);
  
  if(up>histo->GetXaxis()->GetBinLowEdge(histo->GetXaxis()->GetNbins()+1))
    up = histo->GetXaxis()->GetBinLowEdge(histo->GetXaxis()->GetNbins()+1);
  
  const int lowbin = histo->GetXaxis()->FindBin(low);
  const int upbin = histo->GetXaxis()->FindBin(up);
  
  float integral = histo->Integral(lowbin, upbin);
  integral -= histo->GetBinContent(lowbin)*(low-histo->GetXaxis()->GetBinLowEdge(lowbin))/histo->GetXaxis()->GetBinWidth(lowbin);
  integral -= histo->GetBinContent(upbin)*(histo->GetXaxis()->GetBinLowEdge(upbin+1)-up)/histo->GetXaxis()->GetBinWidth(upbin);
  
  return integral;
}