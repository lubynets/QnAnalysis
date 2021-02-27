#include "GraphExtractor.hpp"

#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TH3.h"

std::string StringBinNumber(int number);

int main(int argc, char** argv)
{
  TString filename="/home/user/cbmdir/working/qna/bin_extract/cl.dcmqgsm.apr20.defcuts.nopid.set3.all.root";
  TFile* fileIn = TFile::Open(filename, "read");
  
  GraphExtractor gex;
  Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>* mc_lambda_psi_xx_sgnl = (Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>*)fileIn -> Get("sim/u_sim_PLAIN.Q_psi_PLAIN.x1x1");
  gex.SetDataContainer(mc_lambda_psi_xx_sgnl);
  gex.SetNamesAxesToExclude({"AnaEventHeader_tracks_centrality", "SimParticles_rapidity", "SimParticles_pT"});
  
  std::vector<int> axessizes = gex.GetAxesSizes();
  std::vector<std::vector<double>> axisbinedges = gex.GetAxesBinEdges();
  
  const int C_nbins = axessizes.at(0);
  const int y_nbins = axessizes.at(1);
  const int pT_nbins = axessizes.at(2);
  
  double* C_edges = &axisbinedges.at(0)[0];
  double* y_edges = &axisbinedges.at(1)[0];
  double* pT_edges = &axisbinedges.at(2)[0];
  
  TH3F hv1_mc("hv1_mc", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  hv1_mc.GetXaxis()->SetTitle("centrality, %");
  hv1_mc.GetYaxis()->SetTitle("rapidity");
  hv1_mc.GetZaxis()->SetTitle("p_{T}, GeV");
  
  TFile* fileOut = TFile::Open("out.mcv1.root", "recreate");
  
  for(int iC=0; iC<C_nbins; iC++)
    for(int iy=0; iy<y_nbins; iy++)
      for(int ipT=0; ipT<pT_nbins; ipT++)
      {      
        TGraph* gr = gex.GetGraph({iC, iy, ipT});
        if (gr->GetN() != 1)
          std::cout << "gr->GetN() != 1 !!!\n";
        
        hv1_mc.SetBinContent(iC+1, iy+1, ipT+1, gr->GetPointY(0));
        hv1_mc.SetBinError(iC+1, iy+1, ipT+1, gr->GetErrorY(0));
      } 
  
  hv1_mc.Write();
  
  fileOut -> Close(); 
  
  return 0;
}