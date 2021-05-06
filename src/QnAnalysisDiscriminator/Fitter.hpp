#ifndef Fitter_H
#define Fitter_H

#include "ShapeContainer.hpp"

#include "TH1.h"
#include "TString.h"
#include "TGraph.h"
#include "TF1.h"

class MyFunctor
{
public:  
  
  MyFunctor(ShapeContainer* shape)
  {
    shape_ = shape;
  }
  
  virtual ~MyFunctor() = default;
  
  double signal_fit(double* x, double* par)
  {
    return par[0];
  }
  
  double bckgr_fit(double* x, double* par)
  {
    return par[0] + par[1]*(x[0]-mu);// + par[2]*(x[0]-mu)*(x[0]-mu);
//     return par[0];
    
  }
    
  double operator()(double* x, double* par)
  {
    return (shape_->GetSignal(x[0])*signal_fit(x, par) + shape_->GetBackground(x[0])*bckgr_fit(x, &par[1])) / (shape_->GetSignal(x[0]) + shape_->GetBackground(x[0]));
  }
  
private:
  
  const float mu = 1.11572;   //TODO remove this ugly hardcode ASAP

  ShapeContainer* shape_{nullptr};
};

class Fitter
{
public:
  
  Fitter() = default;
  virtual ~Fitter() = default;  
  
  void SetShape(ShapeContainer* shape) { shape_ = shape; };
  void SetGraphToFit(TGraph* graph) { graph_v_ = graph; };
  double GetVSignal() { return fit_params_.at(0); };
  double GetVSignalError() { return fit_params_errors_.at(0); };
  double GetFitChi2() { return fit_chi2_; };
  int GetFitNdf() { return fit_ndf_; };
  double GetFitChi2Ndf() { return fit_chi2_/fit_ndf_; };
  const std::vector<double>& GetFitParameters() { return fit_params_; };
  const std::vector<double>& GetFitErrors() { return fit_params_errors_; };
  
  void Fit();
  
  
private:
  
  ShapeContainer* shape_{nullptr};
  TGraph* graph_v_{nullptr};
  
  std::vector<double> fit_params_;
  std::vector<double> fit_params_errors_;
  double fit_chi2_{-999.};
  int fit_ndf_{-999};
};

#endif//Fitter_H