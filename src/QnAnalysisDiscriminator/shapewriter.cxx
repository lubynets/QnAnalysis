#include "TString.h"
#include "TFile.h"

#include "ShapeContainer.hpp"

int main(int argc, char** argv)
{
  TString sgnlfilename="/home/user/cbmdir/working/massfit/pfsqa.sgnl_12.root";
  TString bckgrfilename="/home/user/cbmdir/working/massfit/pfsqa.bckgr.root";
  TString v1filename="/home/user/cbmdir/working/massfit/v1graph.toymc.fmd.all.root";
  
  TFile* sgnlfile = TFile::Open(sgnlfilename, "read");
  TH1F* histosgnl = (TH1F*)sgnlfile -> Get("LambdaCandidates/LambdaCandidates_mass");
  TFile* bckgrfile = TFile::Open(bckgrfilename, "read");
  TH1F* histobckgr = (TH1F*)bckgrfile -> Get("LambdaCandidates/LambdaCandidates_mass");

  ShapeContainer sc;
  sc.SetShape(histosgnl, histobckgr);
//   sc.IsA()->SetName("ShapeContainer");
//   sc.IsA()->SetTitle("sctitle");

  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  sc.Write();
  fileOut -> Close();
  
  return 0;
}