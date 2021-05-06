#include "ShapeContainerTensor.hpp"
#include "ShapeFitter.hpp"

#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH3F.h"

#include <string>

std::string StringBinNumber(int number);
float GetChi2TH1FTF1(TH1F* histo, TF1* func);
float GetChi2TH1FTH1F(TH1F* histo1, TH1F* histo2);
void SetAxesNames(TH3F* histo,
                  TString xaxisname="centrality, %",
                  TString yaxisname="rapidity",
                  TString zaxisname="p_{T}, GeV");

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
  
  //***** definition of terms *************************
  //
  // all              shape of sgnl&bckgr together (the only distribution we will operate with in real life)
  // sgnl_mc          shape of sgnl (MC-true)
  // bckgr_mc         shape of bckgr (MC-bckgr)
  // sgnl_mc_fit      fit of sgnl_mc
  // bckgr_rec        exclude peak region from all
  // bckgr_rec_fit    fit of bckgr_rec
  // sgnl_rec         all minus bckgr_fit
  // sgnl_rec_fit     fit of sgnl_rec
  // We consider 2 types of chi2's: chi2 of fit and chi2 of difference between 2 histograms / functions / histogram and function
  //
  //***************************************************
  
  TH3F hchi2_bckgr_rec_fit("hchi2_bckgr_rec_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  SetAxesNames(&hchi2_bckgr_rec_fit);
  
  TH3F hchi2_bckgr_mc_and_bckgr_rec_fit("hchi2_bckgr_mc_and_bckgr_rec_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  SetAxesNames(&hchi2_bckgr_mc_and_bckgr_rec_fit);
  
  TH3F hchi2_sgnl_mc_and_sgnl_rec("hchi2_sgnl_mc_and_sgnl_rec", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  SetAxesNames(&hchi2_sgnl_mc_and_sgnl_rec);
    
  TH3F hchi2_sgnl_mc_fit("hchi2_sgnl_mc_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  SetAxesNames(&hchi2_sgnl_mc_fit);
  
  TH3F hchi2_sgnl_rec_fit("hchi2_sgnl_rec_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  SetAxesNames(&hchi2_sgnl_rec_fit);
  
  TH3F hchi2_sgnl_mc_and_sgnl_rec_fit("hchi2_sgnl_mc_and_sgnl_rec_fit", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  SetAxesNames(&hchi2_sgnl_mc_and_sgnl_rec_fit);
  
  ShapeContainerTensor sct;
  sct.SetFrame({C_nbins, y_nbins, pT_nbins});
  
  TFile* fileOut = TFile::Open("shapetensor_fit.root", "recreate");
  fileOut -> mkdir("bckgr_mc_and_bckgr_rec_fit"); // is chi2 between
  fileOut -> mkdir("sgnl_mc_and_sgnl_rec");       // is chi2 between
  fileOut -> mkdir("sgnl_mc_and_sgnl_mc_fit");    // is chi2 of fit  
  fileOut -> mkdir("sgnl_rec_and_sgnl_rec_fit");  // is chi2 of fit
  fileOut -> mkdir("sgnl_mc_and_sgnl_rec_fit");   // is chi2 between
  
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
                
        fileOut -> cd("bckgr_mc_and_bckgr_rec_fit");
        TCanvas c1("", "", 1500, 900);
        c1.cd();
        histobckgr -> Draw();
        sftr.GetFuncBckgr() -> Draw("same");
        c1.Write(binname.c_str());
        
        fileOut -> cd("sgnl_mc_and_sgnl_rec");
        TCanvas c2("", "", 1500, 900);
        c2.cd();
        histosgnl -> Draw();
        sftr.GetHistoSgnl() -> SetLineColor(kRed);
        sftr.GetHistoSgnl() -> Draw("same");
        c2.Write(binname.c_str());
        
        fileOut -> cd("sgnl_mc_and_sgnl_mc_fit");
        TCanvas c3("", "", 1500, 900);
        c3.cd();
        histosgnl -> Draw();
        ShapeFitter sftr_sgnl_mc(histosgnl);  // can be initialized with any histogram, not important
        TF1* func_sgnl_fitmc = sftr_sgnl_mc.FitSgnl(histosgnl, ShapeFitter::mu - 15*ShapeFitter::sigma, ShapeFitter::mu + 15*ShapeFitter::sigma);
        func_sgnl_fitmc -> Draw("same");
        c3.Write(binname.c_str());
        
        fileOut -> cd("sgnl_rec_and_sgnl_rec_fit");
        TCanvas c4("", "", 1500, 900);
        c4.cd();
        sftr.GetHistoSgnl() -> SetLineColor(kBlue);
        sftr.GetHistoSgnl() -> Draw();
        sftr.GetFuncSgnl() -> Draw("same");
        c4.Write(binname.c_str());
        
        fileOut -> cd("sgnl_mc_and_sgnl_rec_fit");
        TCanvas c5("", "", 1500, 900);
        c5.cd();
        histosgnl -> Draw();
        sftr.GetFuncSgnl() -> Draw("same");
        c5.Write(binname.c_str());
        
        hchi2_bckgr_rec_fit.SetBinContent(iC+1, iy+1, ipT+1, sftr.GetChi2BckgrFit());
        hchi2_bckgr_mc_and_bckgr_rec_fit.SetBinContent(iC+1, iy+1, ipT+1, GetChi2TH1FTF1(histobckgr, sftr.GetFuncBckgr()));
        hchi2_sgnl_mc_and_sgnl_rec.SetBinContent(iC+1, iy+1, ipT+1, GetChi2TH1FTH1F(histosgnl, sftr.GetHistoSgnl()));
        hchi2_sgnl_mc_fit.SetBinContent(iC+1, iy+1, ipT+1, func_sgnl_fitmc->GetChisquare() / func_sgnl_fitmc->GetNDF());
        hchi2_sgnl_rec_fit.SetBinContent(iC+1, iy+1, ipT+1, sftr.GetChi2SgnlFit());
        hchi2_sgnl_mc_and_sgnl_rec_fit.SetBinContent(iC+1, iy+1, ipT+1, GetChi2TH1FTF1(histosgnl, sftr.GetFuncSgnl()));
      }
  
  fileOut -> cd();
  sct.Write("shapetensor");
  hchi2_bckgr_rec_fit.Write();
  hchi2_bckgr_mc_and_bckgr_rec_fit.Write();
  hchi2_sgnl_mc_and_sgnl_rec.Write();
  hchi2_sgnl_mc_fit.Write();
  hchi2_sgnl_rec_fit.Write();
  hchi2_sgnl_mc_and_sgnl_rec_fit.Write();
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
    if(histo->GetBinError(iBin) == 0.) continue;
    const float delta = (func->Eval(histo->GetBinCenter(iBin)) - histo->GetBinContent(iBin)) / histo->GetBinError(iBin);
    chi2 += delta*delta;
    ndf++;
  }
  
  ndf -= func->GetNpar();
  
  std::cout << "chi2/ndf = " << chi2 << " / " << ndf << "\n";
  
  return chi2/ndf;
}

float GetChi2TH1FTH1F(TH1F* histo1, TH1F* histo2)
{
  int ndf = 0;
  float chi2 = 0.f;
  for(int iBin=1; iBin<=histo1->GetNbinsX(); iBin++)
  {
    if(histo1->GetBinError(iBin) == 0. || histo2->GetBinError(iBin) == 0.) continue;
    const float v1 = histo1->GetBinContent(iBin);
    const float v2 = histo2->GetBinContent(iBin);
    const float e1 = histo1->GetBinError(iBin);
    const float e2 = histo2->GetBinError(iBin);
    chi2 += (v1-v2)*(v1-v2)/(e1*e1 + e2*e2);
    ndf++;
  }
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

void SetAxesNames(TH3F* histo, TString xaxisname, TString yaxisname, TString zaxisname)
{
  histo -> GetXaxis() -> SetTitle(xaxisname);
  histo -> GetYaxis() -> SetTitle(yaxisname);
  histo -> GetZaxis() -> SetTitle(zaxisname);
}