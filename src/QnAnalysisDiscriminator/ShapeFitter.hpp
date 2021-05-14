#ifndef ShapeFitter_H
#define ShapeFitter_H

#include "TH1F.h"
#include "TF1.h"

class MyFunctorShape
{
public:
  
  MyFunctorShape() = default;
  virtual ~MyFunctorShape() = default;
  
  double operator()(double* x, double* par);
  
private:
  
  static constexpr float x0_internal_ = 0.9e-3;   // for lambda function need static constexpr instead of const. WHY?
  static constexpr float x0_external_ = 1.3e-3;
};

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
  
  TF1* FitSgnl(TH1F* histo, float left, float right) const;           // TODO move to private after debugging
  
private:
  
  TH1F* ExcludeInterval(TH1F* histo, float left, float right) const;
  TF1* FitBckgr(TH1F* histo, float left, float right) const;
  TH1F* SubtractBckgr(TH1F* histo, TF1* func, float left, float right) const;
    
  TH1F* histo_all_{nullptr};
  TH1F* histo_sgnl_{nullptr};
  TF1* func_bckgr_{nullptr};
  TF1* func_sgnl_{nullptr};
  float chi2_bckgr_fit_{-799.};
  float chi2_sgnl_fit_{-799.};
  
  const float mu = 1.115683;                     //TODO remove this hardcode
  const float sigma = 0.00145786;  
};


#endif // ShapeFitter_H