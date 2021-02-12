#include "NDimFrame.hpp"

#include <iostream>
#include <stdexcept>

int NDimFrame::GetGlobalIndex(std::vector<int> i)
{
  if(n_.size() != i.size())
    throw std::runtime_error("Number of sides is not equal to dimensionality of n-D cube");
  for(int j=0; j<n_.size(); j++)
    if(i.at(j) < 0 || i.at(j)>=n_.at(j))
      throw std::runtime_error("Wrong index number");    
    
  return GetGlobalIndexInternal(i, n_);
}

int NDimFrame::GetGlobalIndexInternal(std::vector<int> i, std::vector<int> n)
{
  if(n.size() == 1)
    return i.back();
  else
  {
    const int nback = n.back();
    const int iback = i.back();
    i.pop_back();
    n.pop_back();
    return nback*NDimFrame::GetGlobalIndexInternal(i, n) + iback;
  }
}