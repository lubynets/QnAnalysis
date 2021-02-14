#ifndef GraphExtractor_H
#define GraphExtractor_H

#include "DataContainer.hpp"

#include "TGraph.h"

#include <vector>
#include <string>

class GraphExtractor
{
public:
  
  GraphExtractor() = default;
  virtual ~GraphExtractor() = default;
  
  void SetDataContainer(Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>* datacontainer) { data_container_ = datacontainer; };
  void SetNamesAxesToExclude(std::vector<std::string> names_axes_to_exclude) { names_axes_to_exclude_ = names_axes_to_exclude; };
  
  TGraph* GetGraph(std::vector<int> bin);
  
private:
  
  Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>* data_container_{nullptr};
  std::vector<std::string> names_axes_to_exclude_;
};

#endif  //GraphExtractor_H