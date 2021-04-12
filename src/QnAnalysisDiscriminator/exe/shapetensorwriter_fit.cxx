#include "ShapeContainerTensor.hpp"
#include "ShapeFitter.hpp"

#include "TString.h"
#include "TFile.h"

#include <string>

std::string StringBinNumber(int number);

int main(int argc, char** argv)
{
  TString allfilename="/home/user/cbmdir/working/qna/shapes/out.mass3D.apr20.dcmqgsm.nopid.lightcuts1.set4.all.root";
  
  TFile* allfile = TFile::Open(allfilename, "read");
  
  TH1F* histoall = nullptr;
  TF1* funcsgnl = nullptr;
  
  const int C_nbins = 3;
  const int y_nbins = 4;
  const int pT_nbins = 4;
  
  ShapeContainerTensor sct;
  sct.SetFrame({C_nbins, y_nbins, pT_nbins});
  
  for(int iC=0; iC<C_nbins; iC++)
    for(int iy=0; iy<y_nbins; iy++)
      for(int ipT=0; ipT<pT_nbins; ipT++)
      {
        std::string binname = "C" + StringBinNumber(iC+1) + "_y" + StringBinNumber(iy+1) + "_pT" + StringBinNumber(ipT+1);
        histoall = (TH1F*) allfile -> Get(binname.c_str());
        ShapeFitter sftr(histoall);
        sftr.Fit();
        sct.SetShape(sftr.GetHistoSgnl(), sftr.GetFuncBckgr(), {iC, iy, ipT});
      }
  
  TFile* fileOut = TFile::Open("shapetensor_fit.root", "recreate");
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