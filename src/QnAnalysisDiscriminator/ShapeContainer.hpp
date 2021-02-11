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

  float GetSignal(float x);
  float GetBackground(float x);  
  
private:
  
  TH1F* histo_sgnl_{nullptr};
  TH1F* histo_bckgr_{nullptr};
  TF1* func_sgnl_{nullptr};
  TF1* func_bckgr_{nullptr};
  
  bool is_histo_{false};        //TODO generalize for combining TH1 & TF1
  
  ClassDef(ShapeContainer, 1);
};

#endif//ShapeContainer_H