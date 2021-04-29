#ifndef ShapeFitter_H
#define ShapeFitter_H

#include "TH1F.h"
#include "TF1.h"

class ShapeFitter
{
public:
  
  ShapeFitter(TH1F* histo);
  virtual ~ShapeFitter() = default;
  
  TH1F* GetHistoSgnl(){ return histo_sgnl_; };
  TF1* GetFuncBckgr(){ return func_bckgr_; };
  void Fit();
  float GetChi2BckgrFit() { return chi2_bckgr_fit_ ; };
  
private:
  
  TH1F* ExcludeInterval(TH1F* histo, float left, float right);
  TF1* FitBckgr(TH1F* histo, float left, float right);
  TH1F* SubtractBckgr(TH1F* histo, TF1* func, float left, float right);
  
  TH1F* histo_all_{nullptr};
  TH1F* histo_sgnl_{nullptr};
  TF1* func_bckgr_{nullptr};  
  float chi2_bckgr_fit_{-799.};
};
#endif // ShapeFitter_H