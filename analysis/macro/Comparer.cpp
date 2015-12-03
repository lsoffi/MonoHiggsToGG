#include "Comparer.hh"

Comparer::Comparer( SamplePairVec Samples, const ColorMap colorMap, const Double_t inLumi, const DblVec puweights, const TString indir, const TString outdir, const Bool_t Blind, const TString type){

  fType = type;
  lumi = inLumi;
  fInDir = indir;
  fOutDir = outdir;
  doBlind = Blind;
  TString fOut = "all"; 
  fPUWeights = puweights;

  // Make output file
  MakeOutDirectory(Form("%s%s",fOutDir.Data(),fOut.Data()));
  fOutFile = new TFile(Form("%s%s/combplots.root",fOutDir.Data(),fOut.Data()),"RECREATE");
  CheckValidFile(fOutFile, Form("%s%s/combplots.root",fOutDir.Data(),fOut.Data())); 

  // Read samples
  for (SamplePairVecIter iter = Samples.begin(); iter != Samples.end(); ++iter){
    if ( (*iter).second == 1 ) {fBkgNames.push_back((*iter).first);}      // background
    else if ( (*iter).second == 0 ) {fSigNames.push_back((*iter).first);} // signal
    else {fDataNames.push_back((*iter).first);}			          // data
  }

  fNData = fDataNames.size();
  fNBkg = fBkgNames.size();
  fNSig = fSigNames.size();

  // Initialize variables for the tree
  Comparer::InitVariables();
  fNTH1D = fTH1DNames.size();

  // Open input files & make the TChain w/ trees
  Comparer::GetInFilesAndMakeTChain();
  Comparer::SetTChainBranchAddresses();
  nentries = fAllChain->GetEntries();
  
  // Check that TChain is filled
  if (nentries == 0){
    std::cout << "nentries in TChain is 0! Exiting." << std::endl;
    return;
  }

  // Define colorMap & title
  fColorMap = colorMap;

  fSampleTitleMap["DoubleEG"]		= "Data";
  fSampleTitleMap["QCD"] 		= "QCD";
  fSampleTitleMap["GJets"]		= "#gamma + Jets";
  fSampleTitleMap["VH"]			= "V + H";
  fSampleTitleMap["DYJetsToLL"]		= "Drell-Yan";
  fSampleTitleMap["GluGluHToGG"]	= "H #rightarrow #gamma#gamma (ggH)";
  fSampleTitleMap["DiPhoton"]		= "#gamma + #gamma";
  fSampleTitleMap["DMHtoGG_M1"]		= "m_{#chi} = 1 GeV";//#bar{#chi}#chi HH ,m_{#chi} = 1 GeV";
  fSampleTitleMap["DMHtoGG_M10"]	= "m_{#chi} = 10 GeV";//#bar{#chi}#chi HH ,m_{#chi} = 10 GeV";
  fSampleTitleMap["DMHtoGG_M100"]	= "m_{#chi} = 100 GeV";//#bar{#chi}#chi HH ,m_{#chi} = 100 GeV";
  fSampleTitleMap["DMHtoGG_M1000"]	= "m_{#chi} = 1000 GeV";//#bar{#chi}#chi HH ,m_{#chi} = 1000 GeV";
  fSampleTitleMap["FakeData"]		= "FakeData";  
  fSampleTitleMap["FakeDataII"]		= "Test";
  fSampleTitleMap["2HDM_mZP600"]	= "m_{Z'} = 600 GeV";//#bar{#chi}#chi HH ,m_{#chi} = 1 GeV";
  fSampleTitleMap["2HDM_mZP800"]	= "m_{Z'} = 800 GeV";//#bar{#chi}#chi HH ,m_{#chi} = 1 GeV";
  fSampleTitleMap["2HDM_mZP1000"]	= "m_{Z'} = 1000 GeV";//#bar{#chi}#chi HH ,m_{#chi} = 1 GeV";
  fSampleTitleMap["2HDM_mZP1200"]	= "m_{Z'} = 1200 GeV";//#bar{#chi}#chi HH ,m_{#chi} = 10 GeV";
  fSampleTitleMap["2HDM_mZP1400"]	= "m_{Z'} = 1400 GeV";//#bar{#chi}#chi HH ,m_{#chi} = 10 GeV";
  fSampleTitleMap["2HDM_mZP1700"]	= "m_{Z'} = 1700 GeV";//#bar{#chi}#chi HH ,m_{#chi} = 100 GeV";
  fSampleTitleMap["2HDM_mZP2000"]	= "m_{Z'} = 2000 GeV";//#bar{#chi}#chi HH ,m_{#chi} = 1000 GeV";
  fSampleTitleMap["2HDM_mZP2500"]	= "m_{Z'} = 2500 GeV";//#bar{#chi}#chi HH ,m_{#chi} = 1000 GeV";


}// end Comparer::Comparer

Comparer::~Comparer(){
  std::cout << "Finished & Deleting" << std::endl;

  for (UInt_t data = 0; data < fNData; data++) { delete fDataFiles[data]; }
  for (UInt_t mc = 0; mc < fNBkg; mc++) { delete fBkgFiles[mc]; }
  for (UInt_t mc = 0; mc < fNSig; mc++) { delete fSigFiles[mc]; }

  delete fOutFile;
  
}// end Comparer::~Comparer

void Comparer::DoComparison(){
 
  Double_t wgt = 1;
 
  for (UInt_t entry = 0; entry < nentries; entry++){
    fAllChain->GetEntry(entry);

    //get weight for each event
    if (sampleID > 0 && sampleID < 10000){// MC
      wgt = (weight)*fPUWeights[nvtx];
    }
    else wgt = 1.;// don't apply weight to data

    //apply selection
    if (hltDiphoton30Mass95 == 1){
      if (mgg >= 100 && mgg <= 200){
         if (sampleID >= 10000 && doBlind){
           
         }
         else{

         }
 
      }// end mgg selection
    }// end trigger selection 

  }// end loop over entries in TChain

}// end Comparer::DoComparison

void Comparer::GetInFilesAndMakeTChain(){
 // open input files into TFileVec for data
 fDataFiles.resize(fNData);
 TStrVec datafile;
 datafile.resize(fNData);
 for (UInt_t data = 0; data < fNData; data++){
   datafile[data] = Form("%s/%s.root",fInDir.Data(),fDataNames[data].Data());
   fDataFiles[data] = TFile::Open(datafile[data].Data());
   CheckValidFile(fDataFiles[data],datafile[data]);
 }

 // open input files into TFileVec for bkg
 fBkgFiles.resize(fNBkg);
 TStrVec bkgfile;
 bkgfile.resize(fNBkg);
 for (UInt_t mc = 0; mc < fNBkg; mc++) {
   bkgfile[mc] = Form("%s/%s.root",fInDir.Data(),fBkgNames[mc].Data());
   fBkgFiles[mc] = TFile::Open(bkgfile[mc].Data());
   CheckValidFile(fBkgFiles[mc],bkgfile[mc]);
 }

 // open input files into TFileVec for sig
 fSigFiles.resize(fNSig);
 TStrVec sigfile;
 sigfile.resize(fNSig);
 for (UInt_t mc = 0; mc < fNSig; mc++) {
   sigfile[mc] = Form("%s/%s.root",fInDir.Data(),fSigNames[mc].Data());
   fSigFiles[mc] = TFile::Open(sigfile[mc].Data());
   CheckValidFile(fSigFiles[mc],sigfile[mc]);
 }

 // make a TChain that has all of the trees from all samples
 // but first check that all of the DiPhotonTree trees are found
 // can get information from each sample by looking at the sampleID
 fAllChain = new TChain("DiPhotonTree");
 for (UInt_t data = 0; data < fNData; data++){
   TTree * tpho = (TTree*)fDataFiles[data]->Get("DiPhotonTree"); 
   CheckValidTree(tpho,"DiPhotonTree",Form("%s/%s.root",fInDir.Data(),fDataNames[data].Data())); 
   fAllChain->Add(datafile[data]);
 }
 for (UInt_t mc = 0; mc < fNBkg; mc++){
   TTree * tpho = (TTree*)fBkgFiles[mc]->Get("DiPhotonTree"); 
   CheckValidTree(tpho,"DiPhotonTree",Form("%s/%s.root",fInDir.Data(),fBkgNames[mc].Data())); 
   fAllChain->Add(bkgfile[mc]);
 }
 for (UInt_t mc = 0; mc < fNSig; mc++){
   TTree * tpho = (TTree*)fSigFiles[mc]->Get("DiPhotonTree"); 
   CheckValidTree(tpho,"DiPhotonTree",Form("%s/%s.root",fInDir.Data(),fSigNames[mc].Data())); 
   fAllChain->Add(sigfile[mc]);
 }

}// end Comparer::GetInFilesAndMakeTChain


TH1D * Comparer::MakeTH1DPlot(const TString hname, const TString htitle, const Int_t nbins, const Double_t xlow, const Double_t xhigh, const TString xtitle, const TString ytitle){
  TH1D * hist = new TH1D(hname.Data(),htitle.Data(),nbins,xlow,xhigh);
  hist->GetXaxis()->SetTitle(xtitle.Data());
  hist->GetYaxis()->SetTitle(ytitle.Data());
  hist->Sumw2();
  gStyle->SetOptStat(1111111);
  return hist;
}// end Comparer::MakeTH1DPlot



void Comparer::InitVariables(){
  fTH1DNames.push_back("mgg");


}// end Comparer::InitVariables

void Comparer::SetTChainBranchAddresses(){

   fAllChain->SetBranchAddress("xsecWeight", &xsecWeight, &b_xsecWeight);
   fAllChain->SetBranchAddress("weight", &weight, &b_weight);
   fAllChain->SetBranchAddress("mggNominal", &mggNominal, &b_mggNominal);
   fAllChain->SetBranchAddress("mggGen", &mggGen, &b_mggGen);
   fAllChain->SetBranchAddress("run", &run, &b_run);
   fAllChain->SetBranchAddress("event", &event, &b_event);
   fAllChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fAllChain->SetBranchAddress("rho", &rho, &b_rho);
   fAllChain->SetBranchAddress("sampleID", &sampleID, &b_sampleID);
   fAllChain->SetBranchAddress("totXsec", &totXsec, &b_totXsec);
   fAllChain->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
   fAllChain->SetBranchAddress("pu_n", &pu_n, &b_pu_n);
   fAllChain->SetBranchAddress("sumDataset", &sumDataset, &b_sumDataset);
   fAllChain->SetBranchAddress("perEveW", &perEveW, &b_perEveW);
   fAllChain->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
   fAllChain->SetBranchAddress("pfmetPhi", &pfmetPhi, &b_pfmetPhi);
   fAllChain->SetBranchAddress("pfmetSumEt", &pfmetSumEt, &b_pfmetSumEt);
   fAllChain->SetBranchAddress("t1pfmet", &t1pfmet, &b_t1pfmet);
   fAllChain->SetBranchAddress("t1pfmetPhi", &t1pfmetPhi, &b_t1pfmetPhi);
   fAllChain->SetBranchAddress("t1pfmetSumEt", &t1pfmetSumEt, &b_t1pfmetSumEt);
   fAllChain->SetBranchAddress("calomet", &calomet, &b_calomet);
   fAllChain->SetBranchAddress("calometPhi", &calometPhi, &b_calometPhi);
   fAllChain->SetBranchAddress("calometSumEt", &calometSumEt, &b_calometSumEt);
   fAllChain->SetBranchAddress("ptgg", &ptgg, &b_ptgg);
   fAllChain->SetBranchAddress("mgg", &mgg, &b_mgg);
   fAllChain->SetBranchAddress("eventClass", &eventClass, &b_eventClass);
   fAllChain->SetBranchAddress("pt1", &pt1, &b_pt1);
   fAllChain->SetBranchAddress("ptOverM1", &ptOverM1, &b_ptOverM1);
   fAllChain->SetBranchAddress("eta1", &eta1, &b_eta1);
   fAllChain->SetBranchAddress("phi1", &phi1, &b_phi1);
   fAllChain->SetBranchAddress("sceta1", &sceta1, &b_sceta1);
   fAllChain->SetBranchAddress("r91", &r91, &b_r91);
   fAllChain->SetBranchAddress("sieie1", &sieie1, &b_sieie1);
   fAllChain->SetBranchAddress("hoe1", &hoe1, &b_hoe1);
   fAllChain->SetBranchAddress("scRawEne1", &scRawEne1, &b_scRawEne1);
   fAllChain->SetBranchAddress("chiso1", &chiso1, &b_chiso1);
   fAllChain->SetBranchAddress("phoiso1", &phoiso1, &b_phoiso1);
   fAllChain->SetBranchAddress("neuiso1", &neuiso1, &b_neuiso1);
   fAllChain->SetBranchAddress("pt2", &pt2, &b_pt2);
   fAllChain->SetBranchAddress("ptOverM2", &ptOverM2, &b_ptOverM2);
   fAllChain->SetBranchAddress("eta2", &eta2, &b_eta2);
   fAllChain->SetBranchAddress("phi2", &phi2, &b_phi2);
   fAllChain->SetBranchAddress("sceta2", &sceta2, &b_sceta2);
   fAllChain->SetBranchAddress("r92", &r92, &b_r92);
   fAllChain->SetBranchAddress("sieie2", &sieie2, &b_sieie2);
   fAllChain->SetBranchAddress("hoe2", &hoe2, &b_hoe2);
   fAllChain->SetBranchAddress("scRawEne2", &scRawEne2, &b_scRawEne2);
   fAllChain->SetBranchAddress("chiso2", &chiso2, &b_chiso2);
   fAllChain->SetBranchAddress("phoiso2", &phoiso2, &b_phoiso2);
   fAllChain->SetBranchAddress("neuiso2", &neuiso2, &b_neuiso2);
   fAllChain->SetBranchAddress("eleveto1", &eleveto1, &b_eleveto1);
   fAllChain->SetBranchAddress("eleveto2", &eleveto2, &b_eleveto2);
   fAllChain->SetBranchAddress("presel1", &presel1, &b_presel1);
   fAllChain->SetBranchAddress("presel2", &presel2, &b_presel2);
   fAllChain->SetBranchAddress("sel1", &sel1, &b_sel1);
   fAllChain->SetBranchAddress("sel2", &sel2, &b_sel2);
   fAllChain->SetBranchAddress("tightsel1", &tightsel1, &b_tightsel1);
   fAllChain->SetBranchAddress("tightsel2", &tightsel2, &b_tightsel2);
   fAllChain->SetBranchAddress("genmatch1", &genmatch1, &b_genmatch1);
   fAllChain->SetBranchAddress("genmatch2", &genmatch2, &b_genmatch2);
   fAllChain->SetBranchAddress("geniso1", &geniso1, &b_geniso1);
   fAllChain->SetBranchAddress("geniso2", &geniso2, &b_geniso2);
   fAllChain->SetBranchAddress("vtxIndex", &vtxIndex, &b_vtxIndex);
   fAllChain->SetBranchAddress("vtxX", &vtxX, &b_vtxX);
   fAllChain->SetBranchAddress("vtxY", &vtxY, &b_vtxY);
   fAllChain->SetBranchAddress("vtxZ", &vtxZ, &b_vtxZ);
   fAllChain->SetBranchAddress("genVtxX", &genVtxX, &b_genVtxX);
   fAllChain->SetBranchAddress("genVtxY", &genVtxY, &b_genVtxY);
   fAllChain->SetBranchAddress("genVtxZ", &genVtxZ, &b_genVtxZ);
   fAllChain->SetBranchAddress("passCHiso1", &passCHiso1, &b_passCHiso1);
   fAllChain->SetBranchAddress("passCHiso2", &passCHiso2, &b_passCHiso2);
   fAllChain->SetBranchAddress("passNHiso1", &passNHiso1, &b_passNHiso1);
   fAllChain->SetBranchAddress("passNHiso2", &passNHiso2, &b_passNHiso2);
   fAllChain->SetBranchAddress("passPHiso1", &passPHiso1, &b_passPHiso1);
   fAllChain->SetBranchAddress("passPHiso2", &passPHiso2, &b_passPHiso2);
   fAllChain->SetBranchAddress("passSieie1", &passSieie1, &b_passSieie1);
   fAllChain->SetBranchAddress("passSieie2", &passSieie2, &b_passSieie2);
   fAllChain->SetBranchAddress("passHoe1", &passHoe1, &b_passHoe1);
   fAllChain->SetBranchAddress("passHoe2", &passHoe2, &b_passHoe2);
   fAllChain->SetBranchAddress("passTightCHiso1", &passTightCHiso1, &b_passTightCHiso1);
   fAllChain->SetBranchAddress("passTightCHiso2", &passTightCHiso2, &b_passTightCHiso2);
   fAllChain->SetBranchAddress("passTightNHiso1", &passTightNHiso1, &b_passTightNHiso1);
   fAllChain->SetBranchAddress("passTightNHiso2", &passTightNHiso2, &b_passTightNHiso2);
   fAllChain->SetBranchAddress("passTightPHiso1", &passTightPHiso1, &b_passTightPHiso1);
   fAllChain->SetBranchAddress("passTightPHiso2", &passTightPHiso2, &b_passTightPHiso2);
   fAllChain->SetBranchAddress("passTightSieie1", &passTightSieie1, &b_passTightSieie1);
   fAllChain->SetBranchAddress("passTightSieie2", &passTightSieie2, &b_passTightSieie2);
   fAllChain->SetBranchAddress("passTightHoe1", &passTightHoe1, &b_passTightHoe1);
   fAllChain->SetBranchAddress("passTightHoe2", &passTightHoe2, &b_passTightHoe2);
   fAllChain->SetBranchAddress("hltPhoton26Photon16Mass60", &hltPhoton26Photon16Mass60, &b_hltPhoton26Photon16Mass60);
   fAllChain->SetBranchAddress("hltPhoton36Photon22Mass15", &hltPhoton36Photon22Mass15, &b_hltPhoton36Photon22Mass15);
   fAllChain->SetBranchAddress("hltPhoton42Photon25Mass15", &hltPhoton42Photon25Mass15, &b_hltPhoton42Photon25Mass15);
   fAllChain->SetBranchAddress("hltDiphoton30Mass95", &hltDiphoton30Mass95, &b_hltDiphoton30Mass95);
   fAllChain->SetBranchAddress("hltDiphoton30Mass70", &hltDiphoton30Mass70, &b_hltDiphoton30Mass70);
   fAllChain->SetBranchAddress("hltDiphoton30Mass55", &hltDiphoton30Mass55, &b_hltDiphoton30Mass55);
   fAllChain->SetBranchAddress("hltDiphoton30Mass55PV", &hltDiphoton30Mass55PV, &b_hltDiphoton30Mass55PV);
   fAllChain->SetBranchAddress("hltDiphoton30Mass55EB", &hltDiphoton30Mass55EB, &b_hltDiphoton30Mass55EB);

};// end Comparer::SetTChainBranchAddresses

