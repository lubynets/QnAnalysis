#ifndef NDimFrame_H
#define NDimFrame_H

#include <vector>
#include <stdexcept>

class NDimFrame
{
public:
  
  NDimFrame(std::vector<int> n)         //TODO possibly initialize with Qn::Axes (?)
  {
    n_ = n;
  };
  
  virtual ~NDimFrame() = default;
  
  int GetGlobalIndex(std::vector<int> i);
  
private:
  
  int GetGlobalIndexInternal(std::vector<int> i, std::vector<int> n);
  
  std::vector<int> n_; //vector of sizes (lengths) of cube edges (side)
  
  //TODO possibly add Qn::Axes? OR at least edges?
  
};

#endif//NDimFrame_H