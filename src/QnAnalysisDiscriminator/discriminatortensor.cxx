#include "GraphExtractor.hpp"
#include "ShapeContainerTensor.hpp"
#include "Fitter.hpp"

#include "TFile.h"
#include "TDirectory.h"
#include "TH3.h"

#include <iostream>

std::string StringBinNumber(int number);

int main(int argc, char** argv)
{
  TString shapefilename="/home/user/cbmdir/working/massfit/shapetensor.apr20.dcmqgsm.nopid.defcuts.set3.root";
  TFile* shapefile = TFile::Open(shapefilename, "read");
  ShapeContainerTensor* shcntr = (ShapeContainerTensor*)shapefile -> Get("shapetensor");
  
  TString v1filename="/home/user/cbmdir/working/qna/bin_extract/cl.dcmqgsm.apr20.defcuts.nopid.set3.all.root";
  TFile* v1file = TFile::Open(v1filename, "read");
  
  Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>* lambda_psi_xx = (Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>*)v1file -> Get("rec/RESCALED/u_rec_RESCALED.Q_psi_PLAIN.x1x1");
  
  GraphExtractor gex;
  gex.SetDataContainer(lambda_psi_xx);
  gex.SetNamesAxesToExclude({"AnaEventHeader_tracks_centrality", "RecParticlesMcPid_rapidity", "RecParticlesMcPid_pT"});
         
  std::vector<int> axessizes = gex.GetAxesSizes();
  std::vector<std::vector<double>> axisbinedges = gex.GetAxesBinEdges();
  
  const int C_nbins = axessizes.at(0);
  const int y_nbins = axessizes.at(1);
  const int pT_nbins = axessizes.at(2);
  
  double* C_edges = &axisbinedges.at(0)[0];
  double* y_edges = &axisbinedges.at(1)[0];
  double* pT_edges = &axisbinedges.at(2)[0];
  
  TH3F hsignal("hsignal", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  TH3F hbckgr_0("hbckgr_0", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  TH3F hbckgr_1("hbckgr_1", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  
  TH3F hentries_sgnl("hentries_sgnl", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  TH3F hentries_bckgr("hentries_bckgr", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  
  hsignal.GetXaxis()->SetTitle("centrality, %");
  hsignal.GetYaxis()->SetTitle("rapidity");
  hsignal.GetZaxis()->SetTitle("p_{T}, GeV");
  hbckgr_0.GetXaxis()->SetTitle("centrality, %");
  hbckgr_0.GetYaxis()->SetTitle("rapidity");
  hbckgr_0.GetZaxis()->SetTitle("p_{T}, GeV");
  hbckgr_1.GetXaxis()->SetTitle("centrality, %");
  hbckgr_1.GetYaxis()->SetTitle("rapidity");
  hbckgr_1.GetZaxis()->SetTitle("p_{T}, GeV");
  hentries_sgnl.GetXaxis()->SetTitle("centrality, %");
  hentries_sgnl.GetYaxis()->SetTitle("rapidity");
  hentries_sgnl.GetZaxis()->SetTitle("p_{T}, GeV");
  hentries_bckgr.GetXaxis()->SetTitle("centrality, %");
  hentries_bckgr.GetYaxis()->SetTitle("rapidity");
  hentries_bckgr.GetZaxis()->SetTitle("p_{T}, GeV");
  
  TFile* fileOut = TFile::Open("out.fitter.root", "recreate");
  TDirectory* dirFit = fileOut->mkdir("fit");
  TDirectory* dirPar = fileOut->mkdir("parameters");
  dirFit->cd();
    
  for(int iC=0; iC<C_nbins; iC++)
    for(int iy=0; iy<y_nbins; iy++)
      for(int ipT=0; ipT<pT_nbins; ipT++)
      {
        std::string graphname = "C" + StringBinNumber(iC+1) + "_y" + StringBinNumber(iy+1) + "_pT" + StringBinNumber(ipT+1);
        TGraph* gr = gex.GetGraph({iC, iy, ipT});
        gr -> SetName(graphname.c_str());
        gr -> GetXaxis() -> SetTitle("m_{inv}, GeV");
        
        Fitter fitter;
        fitter.SetShape(shcntr->GetShape({iC, iy, ipT}));
        fitter.SetGraphToFit(gr);
        fitter.Fit();
        
        std::string vsignal = std::to_string(fitter.GetVSignal()) + " pm " + std::to_string(fitter.GetVSignalError());
               
        gr -> SetTitle(vsignal.c_str());
        gr -> Write();
        
        hsignal.SetBinContent(iC+1, iy+1, ipT+1, fitter.GetFitParameters().at(0));
        hsignal.SetBinError(iC+1, iy+1, ipT+1, fitter.GetFitErrors().at(0));
        hbckgr_0.SetBinContent(iC+1, iy+1, ipT+1, fitter.GetFitParameters().at(1));
        hbckgr_0.SetBinError(iC+1, iy+1, ipT+1, fitter.GetFitErrors().at(1));
        hbckgr_1.SetBinContent(iC+1, iy+1, ipT+1, fitter.GetFitParameters().at(2));
        hbckgr_1.SetBinError(iC+1, iy+1, ipT+1, fitter.GetFitErrors().at(2));
        
        hentries_sgnl.SetBinContent(iC+1, iy+1, ipT+1, shcntr->GetShape({iC, iy, ipT})->GetSignalIntegral(1.110, 1.1214));        //TODO remove hardcode
        hentries_bckgr.SetBinContent(iC+1, iy+1, ipT+1, shcntr->GetShape({iC, iy, ipT})->GetBackgroundIntegral(1.110, 1.1214)); 
      }
      
  dirPar->cd();
  hsignal.Write();
  hbckgr_0.Write();
  hbckgr_1.Write();
  
  hentries_sgnl.Write();
  hentries_bckgr.Write();
  
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