#ifndef ShapeFitter_H
#define ShapeFitter_H

#include "TH1F.h"
#include "TF1.h"

class ShapeFitter
{
public:
  
  ShapeFitter(TH1F* histo);
  virtual ~ShapeFitter() = default;
  
  TH1F* GetHistoSgnl() const { return histo_sgnl_; };
  TF1* GetFuncBckgr() const { return func_bckgr_; };
  TF1* GetFuncSgnl() const { return func_sgnl_; };
  void Fit();
  float GetChi2BckgrFit() const { return chi2_bckgr_fit_ ; };
  float GetChi2SgnlFit() const { return chi2_sgnl_fit_; };
  
private:
  
  TH1F* ExcludeInterval(TH1F* histo, float left, float right) const;
  TF1* FitBckgr(TH1F* histo, float left, float right) const;
  TF1* FitSgnl(TH1F* histo, float left, float right) const;
  TH1F* SubtractBckgr(TH1F* histo, TF1* func, float left, float right) const;
    
  TH1F* histo_all_{nullptr};
  TH1F* histo_sgnl_{nullptr};
  TF1* func_bckgr_{nullptr};
  TF1* func_sgnl_{nullptr};
  float chi2_bckgr_fit_{-799.};
  float chi2_sgnl_fit_{-799.};
  
  const float mu = 1.11572;                     //TODO remove this hardcode
  const float sigma = 0.00145786;  
};
#endif // ShapeFitter_H