#ifndef ShapeContainerTensor_H
#define ShapeContainerTensor_H

#include "ShapeContainer.hpp"
#include "NDimFrame.hpp"

#include "TObject.h"

#include <vector>

class ShapeContainerTensor : public TObject
{
public:
  
  ShapeContainerTensor() = default;
  virtual ~ShapeContainerTensor() = default;
  
  void SetFrame(std::vector<int> lengths);
  void SetShape(TH1F* histosgnl, TH1F* histobckgr, std::vector<int> i);
  void SetShape(TF1* funcsgnl, TF1* funcbckgr, std::vector<int> i);
  void SetShape(TH1F* histosgnl, TF1* funcbckgr, std::vector<int> i);
  void SetShape(TF1* funcsgnl, TH1F* histobckgr, std::vector<int> i);
  void SetChi2BckgrFit(float value, std::vector<int> i) { this -> GetShapeContainer(i) -> SetChi2BckgrFit(value); };
  ShapeContainer* GetShapeContainer(std::vector<int> i) { return shape_container_tensor_ -> at(frame_ -> GetGlobalIndex(i)); };
  float GetSignal(float x, std::vector<int> i) { return this -> GetShapeContainer(i) -> GetSignal(x); };
  float GetBackground(float x, std::vector<int> i) { return this -> GetShapeContainer(i) -> GetBackground(x); };
  float GetChi2BckgrFit(std::vector<int> i) { return this -> GetShapeContainer(i) -> GetChi2BckgrFit(); };
  
private:
  
  std::vector<ShapeContainer*>* shape_container_tensor_{nullptr};
  NDimFrame* frame_{nullptr};
  
  ClassDef(ShapeContainerTensor, 1);
};
//TODO Make this class stable (not sensitive to fine-tuning when initializing it) to avoid segfaults
//TODO Possibly add Qn::Axes info
#endif  //ShapeContainerTensor_H