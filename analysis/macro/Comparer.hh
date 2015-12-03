#ifndef _comparertools_
#define _comparertools_

#include "Style.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TROOT.h"
#include "TMath.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "THStack.h"
#include "TLine.h"

#include <iostream>
#include <cmath>

typedef std::vector<TFile*>   TFileVec;
typedef std::vector<TH1D*>    TH1DVec;
typedef std::vector<TH1DVec>  TH1DVecVec;
typedef std::vector<THStack*> THStackVec;
typedef std::vector<TPad*>    TPadVec;
typedef std::vector<TLegend*> TLegVec;
typedef std::vector<TCanvas*> TCanvVec;
typedef std::vector<TLine*>   TLineVec;
typedef std::vector<Double_t> DblVec;

class Comparer{
public:
  Comparer(const SamplePairVec Samples, const ColorMap colorMap, const Double_t inLumi, const DblVec puweights, const TString inname, const TString outname, const TString type);
  void DoComparison();
  void InitVariables();
  ~Comparer();

private:
  Double_t	lumi;
  TString	fType;

  TFileVec	fDataFiles;
  TFileVec	fBkgFiles;
  TFileVec	fSigFiles;

  TStrMap	fSampleTitleMap;
  ColorMap	fColorMap;

  TString	fOutDir;
  TString	fOutName;
  TFile *	fOutFile;

};
#endif
