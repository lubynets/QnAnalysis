#ifndef ShapeFitter_H
#define ShapeFitter_H

#include "TH1F.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"

class MyFunctorShape
{
public:
  
  MyFunctorShape() = default;
  virtual ~MyFunctorShape() = default;
  
  double operator()(double* x, double* par);
  
private:
  
  const float x0_ = 9e-4;
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
  
// private:
  
  TH1F* ExcludeInterval(TH1F* histo, float left, float right) const;
  TF1* FitBckgr(TH1F* histo, float left, float right) const;
  TF1* FitSgnl(TH1F* histo, float left, float right) const;
  TH1F* SubtractBckgr(TH1F* histo, TF1* func, float left, float right) const;
  
  void Subtract(TMatrixD& m, TMatrixD& v, int from, int what) const;
  void SubtractFromLower(TMatrixD& m, TMatrixD& v, int what) const;
  void SubtractFromUpper(TMatrixD& m, TMatrixD& v, int what) const;
  void NullLDCorner(TMatrixD& m, TMatrixD& v) const;
  void NullRUCorner(TMatrixD& m, TMatrixD& v) const;
  void SetToOnes(TMatrixD& m, TMatrixD& v) const;  
    
  TH1F* histo_all_{nullptr};
  TH1F* histo_sgnl_{nullptr};
  TF1* func_bckgr_{nullptr};
  TF1* func_sgnl_{nullptr};
  float chi2_bckgr_fit_{-799.};
  float chi2_sgnl_fit_{-799.};
  
  static constexpr float mu = 1.11572;                     //TODO remove this hardcode
  static constexpr float sigma = 0.00145786;  
};


#endif // ShapeFitter_H