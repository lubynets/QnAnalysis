#include "ShapeContainerTensor.hpp"

#include "TString.h"
#include "TFile.h"

#include <string>

std::string StringBinNumber(int number);

int main(int argc, char** argv)
{
  TString sgnlfilename="/home/user/cbmdir/working/massfit/out.mass3D.apr20.dcmqgsm.nopid.defcuts.set2.sgnl_12.root";
  TString bckgrfilename="/home/user/cbmdir/working/massfit/out.mass3D.apr20.dcmqgsm.nopid.defcuts.set2.bckgr.root";
  
  TFile* sgnlfile = TFile::Open(sgnlfilename, "read");
  TFile* bckgrfile = TFile::Open(bckgrfilename, "read");
  
  TH1F* histosgnl = nullptr;
  TH1F* histobckgr = nullptr;
  
  const int C_nbins = 3;
  const int y_nbins = 5;
  const int pT_nbins = 5;
  
  ShapeContainerTensor sct;
  sct.SetFrame({C_nbins, y_nbins, pT_nbins});
  
  for(int iC=0; iC<C_nbins; iC++)
    for(int iy=0; iy<y_nbins; iy++)
      for(int ipT=0; ipT<pT_nbins; ipT++)
      {
        std::string binname = "C" + StringBinNumber(iC+1) + "_y" + StringBinNumber(iy+1) + "_pT" + StringBinNumber(ipT+1);
        histosgnl = (TH1F*) sgnlfile -> Get(binname.c_str());
        histobckgr = (TH1F*) bckgrfile -> Get(binname.c_str());
        sct.SetShape(histosgnl, histobckgr, {iC, iy, ipT});
      }
  
  TFile* fileOut = TFile::Open("shapetensor.root", "recreate");
  sct.Write("shapetensor");
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