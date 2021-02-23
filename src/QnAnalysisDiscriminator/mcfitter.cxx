#include "GraphExtractor.hpp"

#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TH3.h"

std::string StringBinNumber(int number);

int main(int argc, char** argv)
{
  TString sgnlfilename="/home/user/cbmdir/working/qna/bin_extract/cl.dcmqgsm.apr20.defcuts.nopid.set2.sgnl_12.root";
  TFile* sgnlfile = TFile::Open(sgnlfilename, "read");
    
  TString bckgrfilename="/home/user/cbmdir/working/qna/bin_extract/cl.dcmqgsm.apr20.defcuts.nopid.set2.bckgr.root";
  TFile* bckgrfile = TFile::Open(bckgrfilename, "read");
  
  GraphExtractor sgnlgex;
  Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>* lambda_psi_xx_sgnl = (Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>*)sgnlfile -> Get("rec/RESCALED/u_rec_RESCALED.Q_psi_PLAIN.x1x1");
  sgnlgex.SetDataContainer(lambda_psi_xx_sgnl);
  sgnlgex.SetNamesAxesToExclude({"AnaEventHeader_tracks_centrality", "RecParticlesMcPid_rapidity", "RecParticlesMcPid_pT"});
  
  GraphExtractor bckgrgex;
  Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>* lambda_psi_xx_bckgr = (Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>*)bckgrfile -> Get("rec/RESCALED/u_rec_RESCALED.Q_psi_PLAIN.x1x1");
  bckgrgex.SetDataContainer(lambda_psi_xx_bckgr);
  bckgrgex.SetNamesAxesToExclude({"AnaEventHeader_tracks_centrality", "RecParticlesMcPid_rapidity", "RecParticlesMcPid_pT"});
  
  std::vector<int> axessizes = sgnlgex.GetAxesSizes();
  std::vector<std::vector<double>> axisbinedges = sgnlgex.GetAxesBinEdges();
  
  const int C_nbins = axessizes.at(0);
  const int y_nbins = axessizes.at(1);
  const int pT_nbins = axessizes.at(2);
  
  double* C_edges = &axisbinedges.at(0)[0];
  double* y_edges = &axisbinedges.at(1)[0];
  double* pT_edges = &axisbinedges.at(2)[0];
  
  TH3F hsignal("hsignal", "HSIGNAL", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  TH3F hbckgr_0("hbckgr_0", "HBCKGR_0", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  TH3F hbckgr_1("hbckgr_1", "HBCKGR_1", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  
  TFile* fileOut = TFile::Open("out.mcfitter.root", "recreate");
  TDirectory* dirSignal = fileOut->mkdir("signalfits");
  TDirectory* dirBckgr = fileOut->mkdir("bckgrfits");
  TDirectory* dirPar = fileOut->mkdir("parameters");
  
  for(int iC=0; iC<C_nbins; iC++)
    for(int iy=0; iy<y_nbins; iy++)
      for(int ipT=0; ipT<pT_nbins; ipT++)
      {
        std::string graphname = "C" + StringBinNumber(iC+1) + "_y" + StringBinNumber(iy+1) + "_pT" + StringBinNumber(ipT+1);
        
        TGraph* grsignal = sgnlgex.GetGraph({iC, iy, ipT});
        grsignal -> SetName(graphname.c_str());
        grsignal -> GetXaxis() -> SetTitle("m_{inv}, GeV");
        
        TF1 fit_signal("fit_signal", "[0]", grsignal->GetPointX(0)-0.01, grsignal->GetPointX(grsignal->GetN()-1));
        grsignal -> Fit("fit_signal");
        
        std::string vsignal = std::to_string(fit_signal.GetParameter(0)) + " pm " + std::to_string(fit_signal.GetParError(0));
        hsignal.SetBinContent(iC+1, iy+1, ipT+1, fit_signal.GetParameter(0));
        hsignal.SetBinError(iC+1, iy+1, ipT+1, fit_signal.GetParError(0));
               
        grsignal -> SetTitle(vsignal.c_str());
        dirSignal -> cd();
        grsignal -> Write();
        //------------------------------------------------------------------
        
        TGraph* grbckgr = bckgrgex.GetGraph({iC, iy, ipT});
        grbckgr -> SetName(graphname.c_str());
        grbckgr -> GetXaxis() -> SetTitle("m_{inv}, GeV");
        
        TF1 fit_bckgr("fit_bckgr", "[0]+[1]*(x-1.11572)", grbckgr->GetPointX(0)-0.01, grbckgr->GetPointX(grbckgr->GetN()-1));
        grbckgr -> Fit("fit_bckgr");
        
//         vsignal = std::to_string(fit_bckgr.GetParameter(0)) + " pm " + std::to_string(fit_bckgr.GetParError(0));
        hbckgr_0.SetBinContent(iC+1, iy+1, ipT+1, fit_bckgr.GetParameter(0));
        hbckgr_0.SetBinError(iC+1, iy+1, ipT+1, fit_bckgr.GetParError(0));
        hbckgr_1.SetBinContent(iC+1, iy+1, ipT+1, fit_bckgr.GetParameter(1));
        hbckgr_1.SetBinError(iC+1, iy+1, ipT+1, fit_bckgr.GetParError(1));
               
//         grbckgr -> SetTitle(vsignal.c_str());
        dirBckgr -> cd();
        grbckgr -> Write();
      }
    
  dirPar -> cd();
  hsignal.Write();
  hbckgr_0.Write();
  hbckgr_1.Write();
  
  fileOut -> Close();  
  
  return 0;
}

std::string StringBinNumber(int number)
{
  if(number<10)
    return "0" + std::to_string(number);
  else
    return std::to_string(number);
}