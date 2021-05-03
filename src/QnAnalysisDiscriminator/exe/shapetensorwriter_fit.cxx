#include "ShapeContainerTensor.hpp"
#include "ShapeFitter.hpp"

#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH3F.h"

#include <string>

std::string StringBinNumber(int number);
float GetChi2TH1FTF1(TH1F* histo, TF1* func);

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
  
  TH3F hchi2_bckgr_func_histo("hchi2_bckgr_func_histo", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  hchi2_bckgr_func_histo.GetXaxis()->SetTitle("centrality, %");
  hchi2_bckgr_func_histo.GetYaxis()->SetTitle("rapidity");
  hchi2_bckgr_func_histo.GetZaxis()->SetTitle("p_{T}, GeV");
  
  TH3F hchi2_sgnl_fit("hchi2_sgnl_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  hchi2_sgnl_fit.GetXaxis()->SetTitle("centrality, %");
  hchi2_sgnl_fit.GetYaxis()->SetTitle("rapidity");
  hchi2_sgnl_fit.GetZaxis()->SetTitle("p_{T}, GeV");
  
  ShapeContainerTensor sct;
  sct.SetFrame({C_nbins, y_nbins, pT_nbins});
  
  TFile* fileOut = TFile::Open("shapetensor_fit.root", "recreate");
  fileOut -> mkdir("bckgr_mc_and_fitrec");
  fileOut -> mkdir("sgnl_mc_and_rec");
  fileOut -> mkdir("sgnl_rec_and_fitrec");
  
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
                
        fileOut -> cd("bckgr_mc_and_fitrec");
        TCanvas c_bckgr("", "", 1500, 900);
        c_bckgr.cd();
        histobckgr -> Draw();
        sftr.GetFuncBckgr() -> Draw("same");
        c_bckgr.Write(binname.c_str());
        
        fileOut -> cd("sgnl_mc_and_rec");
        TCanvas c_sgnl_res("", "", 1500, 900);
        c_sgnl_res.cd();
        histosgnl -> Draw();
        sftr.GetHistoSgnl() -> SetLineColor(kRed);
        sftr.GetHistoSgnl() -> Draw("same");
        c_sgnl_res.Write(binname.c_str());
        
        fileOut -> cd("sgnl_rec_and_fitrec");
        TCanvas c_sgnl_fit("", "", 1500, 900);
        c_sgnl_fit.cd();
        sftr.GetHistoSgnl() -> SetLineColor(kBlue);
        sftr.GetHistoSgnl() -> Draw();
        sftr.GetFuncSgnl() -> Draw("same");
        c_sgnl_fit.Write(binname.c_str());        
        
        hchi2_bckgr_fit.SetBinContent(iC+1, iy+1, ipT+1, sftr.GetChi2BckgrFit());
        hchi2_bckgr_func_histo.SetBinContent(iC+1, iy+1, ipT+1, GetChi2TH1FTF1(histobckgr, sftr.GetFuncBckgr()));
        hchi2_sgnl_fit.SetBinContent(iC+1, iy+1, ipT+1, sftr.GetChi2SgnlFit());
      }
  
  fileOut -> cd();
  sct.Write("shapetensor");
  hchi2_bckgr_fit.Write();
  hchi2_bckgr_func_histo.Write();
  hchi2_sgnl_fit.Write();
  fileOut -> Close();
  
  
  allfile -> Close();
  
  return 0;
}

float GetChi2TH1FTF1(TH1F* histo, TF1* func)
{
  int firstbin = histo -> FindBin(func->GetXmin());
  if(histo->GetBinCenter(firstbin) < func->GetXmin())
    firstbin++;
  
  int lastbin = histo -> FindBin(func->GetXmax());
  if(histo->GetBinCenter(lastbin) > func->GetXmax())
    lastbin--;
  
  int ndf = 0;
  float chi2 = 0.f;
  for(int iBin=firstbin; iBin<=lastbin; iBin++)
  {
    const float delta = (func->Eval(histo->GetBinCenter(iBin)) - histo->GetBinContent(iBin)) / histo->GetBinError(iBin);
    chi2 += delta*delta;
    ndf++;
  }
  
  ndf -= func->GetNpar();
  
  std::cout << "chi2/ndf = " << chi2 << " / " << ndf << "\n";
  
  return chi2/ndf;
}

std::string StringBinNumber(int number)
{
  if(number<10)
    return "0" + std::to_string(number);
  else
    return std::to_string(number);
}