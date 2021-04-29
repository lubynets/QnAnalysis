#include "ShapeContainerTensor.hpp"

ClassImp(ShapeContainerTensor);

void ShapeContainerTensor::SetFrame(std::vector<int> lengths)
{
  int nbins = 1;
  for(auto length : lengths)
  {
    if(length<=0)
      throw std::runtime_error("Length of n-D cube side must have positive value");
    nbins *= length;
  }
  
  frame_ = new NDimFrame(lengths);
  
  shape_container_tensor_ = new std::vector<ShapeContainer*>();
  shape_container_tensor_ -> resize(nbins);
  for(auto& shc : *shape_container_tensor_)
    shc = new ShapeContainer();
}

void ShapeContainerTensor::SetShape(TH1F* histosgnl, TH1F* histobckgr, std::vector<int> i)
{
  this -> GetShapeContainer(i) -> SetShape(histosgnl, histobckgr);
}

void ShapeContainerTensor::SetShape(TF1* funcsgnl, TF1* funcbckgr, std::vector<int> i)
{
  this -> GetShapeContainer(i) -> SetShape(funcsgnl, funcbckgr);
}

void ShapeContainerTensor::SetShape(TH1F* histosgnl, TF1* funcbckgr, std::vector<int> i)
{
  this -> GetShapeContainer(i) -> SetShape(histosgnl, funcbckgr);
}

void ShapeContainerTensor::SetShape(TF1* funcsgnl, TH1F* histobckgr, std::vector<int> i)
{
  this -> GetShapeContainer(i) -> SetShape(funcsgnl, histobckgr);
}