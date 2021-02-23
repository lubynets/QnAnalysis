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

std::vector<int> GraphExtractor::GetAxesSizes()
{
  std::vector<int> sizes;
  for(auto& axisname : names_axes_to_exclude_)
    sizes.push_back(data_container_->GetAxis(axisname.c_str()).GetNBins());
  
  return sizes;
}

std::vector<std::vector<double>> GraphExtractor::GetAxesBinEdges()
{
  std::vector<std::vector<double>> v_binedges;
  for(auto& axisname : names_axes_to_exclude_)
  {
    std::vector<double> binedges;
    for(int i=0; i<=data_container_->GetAxis(axisname.c_str()).GetNBins(); i++)
      binedges.push_back(data_container_->GetAxis(axisname.c_str()).GetLowerBinEdge(i));
    v_binedges.push_back(binedges);
  }
  
  return v_binedges;
}