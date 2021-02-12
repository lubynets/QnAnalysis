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