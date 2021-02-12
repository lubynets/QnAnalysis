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
  ShapeContainer* GetShape(std::vector<int> i) { return shape_container_tensor_ -> at(frame_ -> GetGlobalIndex(i)); };
  float GetSignal(float x, std::vector<int> i) { return shape_container_tensor_ -> at(frame_ -> GetGlobalIndex(i)) -> GetSignal(x); };
  float GetBackground(float x, std::vector<int> i) { return shape_container_tensor_ -> at(frame_ -> GetGlobalIndex(i)) -> GetBackground(x); };
  
private:
  
  std::vector<ShapeContainer*>* shape_container_tensor_{nullptr};
  NDimFrame* frame_{nullptr};
  
  ClassDef(ShapeContainerTensor, 1);
};
//TODO Make this class stable (not sensitive to fine-tuning when initializing it) to avoid segfaults
//TODO Possibly add Qn::Axes info
#endif  //ShapeContainerTensor_H