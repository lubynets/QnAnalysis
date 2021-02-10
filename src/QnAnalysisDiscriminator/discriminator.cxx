#include "Fitter.hpp"

int main(int argc, char** argv)
{
  TString sgnlfilename="/home/user/cbmdir/working/massfit/pfsqa.sgnl_12.root";
  TString bckgrfilename="/home/user/cbmdir/working/massfit/pfsqa.bckgr.root";
  TString v1filename="/home/user/cbmdir/working/massfit/v1graph.toymc.fmd.all.root";
  
  TFile* v1file = TFile::Open(v1filename, "read");
  TGraph* v1graph = (TGraph*)v1file -> Get("v1_invmass");
  
  TFile* sgnlfile = TFile::Open(sgnlfilename, "read");
  TH1F* histosgnl = (TH1F*)sgnlfile -> Get("LambdaCandidates/LambdaCandidates_mass");
  TFile* bckgrfile = TFile::Open(bckgrfilename, "read");
  TH1F* histobckgr = (TH1F*)bckgrfile -> Get("LambdaCandidates/LambdaCandidates_mass");
  
  Fitter fitter;
  fitter.SetSignalShape(histosgnl);
  fitter.SetBackgroundShape(histobckgr);
  fitter.SetGraphToFit(v1graph);
  
  fitter.Fit();
  
  return 0;
}