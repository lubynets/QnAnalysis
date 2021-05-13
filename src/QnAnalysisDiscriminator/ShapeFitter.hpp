#ifndef ShapeFitter_H
#define ShapeFitter_H

#include "TH1F.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TGraphErrors.h"

class MyFunctorShape
{
public:
  
  MyFunctorShape() = default;
  virtual ~MyFunctorShape() = default;
  
  double operator()(double* x, double* par);
  
private:
  
  const float x0_ = 1.1e-3;
};

class ShapeFitter
{
public:
  
  ShapeFitter(TH1F* histo);
  virtual ~ShapeFitter() = default;
  
  TH1F* GetHistoSgnl() const { return histo_sgnl_; };
  TF1* GetFuncBckgr() const { return func_bckgr_.first; };
  TGraphErrors* GetGraphBckgr() const { return graph_bckgr_; };
  TF1* GetFuncSgnl() const { return func_sgnl_; };
  void Fit();
  float GetChi2BckgrFit() const { return chi2_bckgr_fit_ ; };
  float GetChi2SgnlFit() const { return chi2_sgnl_fit_; };
  
// private:
  
  TH1F* ExcludeInterval(TH1F* histo, float left, float right) const;
  std::pair<TF1*, TMatrixDSym*> FitBckgr(TH1F* histo, float left, float right) const;
  TF1* FitSgnl(TH1F* histo, float left, float right) const;
  TH1F* SubtractBckgr(TH1F* histo, TF1* func, float left, float right) const;
  float EvalError(double* x, std::pair<TF1*, TMatrixDSym*> f_and_cov) const;
  TGraphErrors* FuncWithErrors(std::pair<TF1*, TMatrixDSym*> f_and_cov) const;
    
  TH1F* histo_all_{nullptr};
  TH1F* histo_sgnl_{nullptr};
  std::pair<TF1*, TMatrixDSym*> func_bckgr_{nullptr, nullptr};
  TGraphErrors* graph_bckgr_{nullptr};
  TF1* func_sgnl_{nullptr};
  float chi2_bckgr_fit_{-799.};
  float chi2_sgnl_fit_{-799.};
  
  static constexpr float mu = 1.115683;                     //TODO remove this hardcode
  static constexpr float sigma = 0.00145786;  
};
#endif // ShapeFitter_H