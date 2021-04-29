#ifndef ShapeContainer_H
#define ShapeContainer_H

#include "TH1.h"
#include "TF1.h"
#include "TObject.h"

class ShapeContainer : public TObject
{
public:
  
  ShapeContainer() = default;
  virtual ~ShapeContainer() = default;
  
  void SetShape(TH1F* histosgnl, TH1F* histobckgr);
  void SetShape(TF1* funcsgnl, TF1* funcbckgr);
  void SetShape(TH1F* histosgnl, TF1* funcbckgr);
  void SetShape(TF1* funcsgnl, TH1F* histobckgr);
  void SetChi2BckgrFit(float value) { chi2_bckgr_fit_ = value ; };

  float GetSignal(float x);
  float GetBackground(float x);
  float GetChi2BckgrFit() { return chi2_bckgr_fit_ ; };
  
  float GetSignalIntegral(float left, float right);
  float GetBackgroundIntegral(float left, float right);
  
private:
  
  TH1F* histo_sgnl_{nullptr};
  TH1F* histo_bckgr_{nullptr};
  TF1* func_sgnl_{nullptr};
  TF1* func_bckgr_{nullptr};
  float chi2_bckgr_fit_{-999.};
  
  bool is_histo_sgnl_{false};
  bool is_histo_bckgr_{false};
  
  float HistoIntegral(TH1F* histo, float low, float up);
  
  ClassDef(ShapeContainer, 1);
};

#endif//ShapeContainer_H