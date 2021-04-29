#include "ShapeContainerTensor.hpp"
#include "ShapeFitter.hpp"

#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH3F.h"

#include <string>

std::string StringBinNumber(int number);

int main(int argc, char** argv)
{
  TString allfilename="/home/user/cbmdir/working/qna/shapes/out.mass3D.apr20.dcmqgsm.nopid.lightcuts1.set4.all.root";
  TString sgnlfilename="/home/user/cbmdir/working/qna/shapes/out.mass3D.apr20.dcmqgsm.nopid.lightcuts1.set4.sgnl_12.root";
  TString bckgrfilename="/home/user/cbmdir/working/qna/shapes/out.mass3D.apr20.dcmqgsm.nopid.lightcuts1.set4.bckgr.root";
  
  TFile* allfile = TFile::Open(allfilename, "read");
  TFile* sgnlfile = TFile::Open(sgnlfilename, "read");
  TFile* bckgrfile = TFile::Open(bckgrfilename, "read");
  
  TH1F* histoall = nullptr;
  TH1F* histosgnl = nullptr;
  TH1F* histobckgr = nullptr;
  
  const int C_nbins = 3;
  const int y_nbins = 4;
  const int pT_nbins = 4;
  
  double C_edges_array[] = {0, 20, 40, 100};
  double y_edges_array[] = {1.02179, 1.42179, 1.82179, 2.22179, 2.62179};
  double pT_edges_array[] = {0.2, 0.5, 0.8, 1.1, 1.4};
  
  double* C_edges = &C_edges_array[0];
  double* y_edges = &y_edges_array[0];
  double* pT_edges = &pT_edges_array[0];
  
  TH3F hchi2_bckgr_fit("hchi2_bckgr_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  hchi2_bckgr_fit.GetXaxis()->SetTitle("centrality, %");
  hchi2_bckgr_fit.GetYaxis()->SetTitle("rapidity");
  hchi2_bckgr_fit.GetZaxis()->SetTitle("p_{T}, GeV");
  
  ShapeContainerTensor sct;
  sct.SetFrame({C_nbins, y_nbins, pT_nbins});
  
  TFile* fileOut = TFile::Open("shapetensor_fit.root", "recreate");
  fileOut -> mkdir("bckgr");
  fileOut -> mkdir("sgnl");
  
  for(int iC=0; iC<C_nbins; iC++)
    for(int iy=0; iy<y_nbins; iy++)
      for(int ipT=0; ipT<pT_nbins; ipT++)
      {
        std::string binname = "C" + StringBinNumber(iC+1) + "_y" + StringBinNumber(iy+1) + "_pT" + StringBinNumber(ipT+1);
        histoall = (TH1F*) allfile -> Get(binname.c_str());
        histosgnl = (TH1F*) sgnlfile -> Get(binname.c_str());
        histobckgr = (TH1F*) bckgrfile -> Get(binname.c_str());
        ShapeFitter sftr(histoall);
        sftr.Fit();
        sct.SetShape(sftr.GetHistoSgnl(), sftr.GetFuncBckgr(), {iC, iy, ipT});
        sct.SetChi2BckgrFit(sftr.GetChi2BckgrFit(), {iC, iy, ipT});
        
        fileOut -> cd("bckgr");
        TCanvas c_bckgr("", "", 1500, 900);
        c_bckgr.cd();
        histobckgr -> Draw();
        sftr.GetFuncBckgr() -> Draw("same");
        c_bckgr.Write(binname.c_str());
        
        fileOut -> cd("sgnl");
        TCanvas c_sgnl("", "", 1500, 900);
        c_sgnl.cd();
        histosgnl -> Draw();
        sftr.GetHistoSgnl() -> SetLineColor(kRed);
        sftr.GetHistoSgnl() -> Draw("same");
        c_sgnl.Write(binname.c_str());
        
        hchi2_bckgr_fit.SetBinContent(iC+1, iy+1, ipT+1, sftr.GetChi2BckgrFit());
      }
  
  fileOut -> cd();
  sct.Write("shapetensor");
  hchi2_bckgr_fit.Write();
  fileOut -> Close();
  
  
  allfile -> Close();
  
  return 0;
}

std::string StringBinNumber(int number)
{
  if(number<10)
    return "0" + std::to_string(number);
  else
    return std::to_string(number);
}