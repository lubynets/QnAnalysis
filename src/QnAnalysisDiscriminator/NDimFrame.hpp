#ifndef NDimFrame_H
#define NDimFrame_H

#include <vector>

class NDimFrame
{
public:
  
  NDimFrame() = default;
  
  NDimFrame(std::vector<int> n)
  {
    n_ = n;
  };
  
  virtual ~NDimFrame() = default;
  
  int GetGlobalIndex(std::vector<int> i);
  
private:
  
  int GetGlobalIndexInternal(std::vector<int> i, std::vector<int> n);
  
  std::vector<int> n_; //vector of sizes (lengths) of cube edges (side)
};

#endif//NDimFrame_H