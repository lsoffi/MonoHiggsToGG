#include "Comparer.hh"

Comparer::Comparer( SamplePairVec Samples, const ColorMap colorMap, const Double_t inLumi, const DblVec puweights, const TString indir, const TString outdir, const TString type){

  fType = type;
  lumi = inLumi;
  fOutDir = outdir;

  Comparer::InitVariables();

}// end Comparer::Comparer

void Comparer::DoComparison(){

}// end Comparer::DoComparison

Comparer::~Comparer(){
  std::cout << "Finished & Deleting" << std::endl;
}


void Comparer::InitVariables(){

}// end Comparer::InitVariables
