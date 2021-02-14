#include "GraphExtractor.hpp"
#include "ShapeContainerTensor.hpp"
#include "Fitter.hpp"

#include "TFile.h"

std::string StringBinNumber(int number);

int main(int argc, char** argv)
{
  TString shapefilename="/home/user/cbmdir/working/massfit/shapetensor.root";
  TFile* shapefile = TFile::Open(shapefilename, "read");
  ShapeContainerTensor* shcntr = (ShapeContainerTensor*)shapefile -> Get("shapetensor");
  
  TString v1filename="/home/user/cbmdir/working/qna/bin_extract/cl.dcmqgsm.apr20.defcuts.nopid.all.large.root";
  TFile* v1file = TFile::Open(v1filename, "read");
  
  Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>* lambda_psi_xx = (Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>*)v1file -> Get("rec/RESCALED/u_rec_RESCALED.Q_psi_PLAIN.x1x1");
  
  GraphExtractor gex;
  gex.SetDataContainer(lambda_psi_xx);
  gex.SetNamesAxesToExclude({"AnaEventHeader_tracks_centrality", "RecParticlesMcPid_rapidity", "RecParticlesMcPid_pT"});
   
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  
  const int C_nbins = 7;
  const int y_nbins = 7;
  const int pT_nbins = 7;
  
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
      }
  
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