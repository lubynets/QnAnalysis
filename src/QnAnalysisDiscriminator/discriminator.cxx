#include "Fitter.hpp"
#include "TFile.h"

int main(int argc, char** argv)
{
  
  TString shapefilename="/home/user/cbmdir/working/massfit/shape.root";
  TString v1filename="/home/user/cbmdir/working/massfit/v1graph.toymc.fmd.all.root";

  TFile* shapefile = TFile::Open(shapefilename, "read");
  ShapeContainer* shcntr = (ShapeContainer*)shapefile -> Get("shape");
  
  TFile* v1file = TFile::Open(v1filename, "read");
  TGraph* v1graph = (TGraph*)v1file -> Get("v1_invmass");
  
  Fitter fitter;
  fitter.SetShape(shcntr);
  fitter.SetGraphToFit(v1graph);
  
  fitter.Fit();
  
  return 0;
}