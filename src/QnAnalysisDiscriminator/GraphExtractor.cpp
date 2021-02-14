#include "GraphExtractor.hpp"

TGraph* GraphExtractor::GetGraph(std::vector<int> bin)
{
  if(bin.size() != names_axes_to_exclude_.size())
    throw std::runtime_error("Bin vector's dimensionality must equal to number of axes to exclude");
    
  std::vector<Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>> data_container_reduced;
  data_container_reduced.push_back(*data_container_);
  
  for(int j=0; j<bin.size(); j++)
  {
    double left_edge = data_container_ -> GetAxis(names_axes_to_exclude_.at(j).c_str()).GetLowerBinEdge(bin.at(j));
    double right_edge = data_container_ -> GetAxis(names_axes_to_exclude_.at(j).c_str()).GetLowerBinEdge(bin.at(j)+1);
    Qn::Axis<double> axis(names_axes_to_exclude_.at(j), 1, left_edge, right_edge);
    data_container_reduced.push_back(data_container_reduced.back().Select(axis));
  }
  
  TGraph* graph = Qn::DataContainerHelper::ToTGraph(data_container_reduced.back());
  
  return graph;
}