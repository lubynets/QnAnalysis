#ifndef CORRELATION_TASK_H
#define CORRELATION_TASK_H

#include "TChain.h"
#include "TTreeReader.h"

#include "QnTools/ReSampleHelper.hpp"
#include "QnTools/DataContainer.hpp"
#include "QnTools/Axis.hpp"
#include "QnTools/QnDataFrame.hpp"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDataSource.hxx"

#include "QVector.h"
#include "AnalysisSetup.h"

class CorrelationTask {
 public:
   enum Q1Q1_CORRELATION_TYPES{
     SCALAR_PRODUCT=0,
     Q1Q1_EVENT_PLANE,
     u1Q1_EVENT_PLANE
   };
  CorrelationTask() = delete;
  CorrelationTask(const std::string& file, const std::string& treename);
  void AddQ1Q1Correlation(const std::string&a_name, const std::string&b_name, int type=SCALAR_PRODUCT){
    AddQ1Q1Correlations(config_->GetQvectorConfig(a_name), config_->GetQvectorConfig(b_name), type);
  }
  void SetNonZeroOnly(bool non_zero_only) { non_zero_only_ = non_zero_only; }
  void Run();

 private:
  void FillCorrelations();
  void AddQ1Q1Correlations(const Flow::Base::QVector& a, const Flow::Base::QVector& b, int type);
  void AddQ2Q1Q1Correlations(const Flow::Base::QVectorTrack& t, const Flow::Base::QVector& a, const Flow::Base::QVector& b);

  TFile* in_file_{nullptr};
  Flow::Base::AnalysisSetup* config_{nullptr};
  TTree* in_tree_{nullptr};
  TTreeReader reader_;
  ROOT::RDataFrame df_;
  ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager> dfs_;
  Qn::AxisD event_;

  std::map<std::string, ROOT::RDF::RResultPtr<Qn::Correlation::CorrelationActionBase>> correlations_{};

  int verbosity_{2};
  bool non_zero_only_{true};
};

#endif //CORRELATION_TASK_H
