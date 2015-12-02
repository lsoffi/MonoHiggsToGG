#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>

using namespace std;

void fitterFormatting(const char* filename, TString type, TString theSample, TString outFile) {
  cout << "Formatting " << filename << endl;
  cout << "Move to Pasquale's format for fit." << endl;

  TFile *fileOrig = 0;
  TTree *treeOrig = 0;

  fileOrig = TFile::Open(filename);
  if ( fileOrig ) {
    fileOrig->cd();
    TString theTree = TString::Format("DiPhotonTree");
    cout << "Reading " << theTree << endl;
    treeOrig = (TTree*)fileOrig->Get(theTree);
  }
  else {
    cout << "File " << filename << " does not exist!" << endl;
    return;
  }

  fileOrig->cd();
  if (!treeOrig) {
    cout << "Tree does not exist!" << endl;
    return;
  }

  treeOrig->SetMakeClass(0);
  // get number of entries in original tree to loop over
  UInt_t nentriesOrig = treeOrig->GetEntries();

  // vector for new photon categories out
  UInt_t numPhoCat = 5;
  vector<TString> thePhotonCat;
  thePhotonCat.resize(numPhoCat);
  thePhotonCat[0]="all";
  thePhotonCat[1]="EBHighR9";
  thePhotonCat[2]="EBLowR9";
  thePhotonCat[3]="EEHighR9";
  thePhotonCat[4]="EELowR9";

  // vector for new met categories out
  UInt_t numMetCat = 12;
  UInt_t metSpacing = 100; // spacing of met bins in GeV
  vector<TString> theMetCat;
  theMetCat.resize(numMetCat);
  for (UInt_t met=0; met<numMetCat; met++){
     UInt_t metmin = met*metSpacing;
     UInt_t metmax = (met+1)*metSpacing;
     if (met==numMetCat-2) theMetCat[met]=TString::Format("%i",metmin);
     if (met==numMetCat-1) theMetCat[met]="_all";
     else theMetCat[met]=TString::Format("%i_%i",metmin,metmax);
  }

  // make output file and new trees
  cout << "OutputFile: " << outFile << endl;
  TFile *fileNew = TFile::Open(outFile,"RECREATE");
  vector<vector<TTree*> > trees;
  vector<TDirectory*> newDir;
  newDir.resize(numMetCat);  

  // make a tree for each photon category and each met category
  trees.resize(numMetCat);
  for (UInt_t met = 0; met<numMetCat; met++){
    trees[met].resize(numPhoCat);
    TString theNewDir=TString::Format("met%s",theMetCat[met].Data());
    TDirectory *Dir = fileNew->mkdir(theNewDir);
    newDir[met] = Dir;
    newDir[met]->cd();
    for (UInt_t pho = 0; pho<numPhoCat; pho++){
      TString theNewTree = TString::Format("%s_%s_%s",type.Data(),theSample.Data(),thePhotonCat[pho].Data());
      TTree *treeNew = new TTree(theNewTree,theNewTree); 
      trees[met][pho]=treeNew;
    }
  }

  // original tree leaves
  Int_t   run		= 0.;
  Int_t   event		= 0.;   
  Int_t	  nvtx		= 0.;
  Float_t weight	= 0.;
  Float_t rho		= 0.;
  Float_t mgg		= 0.;
  Float_t eta1		= 0.;
  Float_t eta2		= 0.;
  Float_t r91		= 0.;
  Float_t r92		= 0.;
  Float_t t1pfmet	= 0.;
  Float_t t1pfmetPhi	= 0.;
  Float_t pt1		= 0.;
  Float_t pt2		= 0.;
  Float_t ptgg		= 0.;
  Float_t chiso1 	= 0.;
  Float_t chiso2 	= 0.;
  Float_t phoiso1	= 0.;
  Float_t phoiso2	= 0.;
  Float_t neuiso1 	= 0.;
  Float_t neuiso2 	= 0.;
  Float_t phi1		= 0.;  
  Float_t phi2		= 0.;  
  Float_t hoe1		= 0.;
  Float_t hoe2		= 0.;
  Float_t sieie1 	= 0.;
  Float_t sieie2 	= 0.;
  Float_t sceta1	= 0.;
  Float_t sceta2	= 0.;
  Int_t   eleveto1	= 0.;
  Int_t   eleveto2	= 0.;

  // branches from original tree
  TBranch *b_run;
  TBranch *b_event;
  TBranch *b_nvtx;
  TBranch *b_rho;
  TBranch *b_weight;
  TBranch *b_mgg;
  TBranch *b_ptgg;
  TBranch *b_eta1;
  TBranch *b_eta2;
  TBranch *b_r91;
  TBranch *b_r92;
  TBranch *b_pt1;
  TBranch *b_pt2;
  TBranch *b_t1pfmet;
  TBranch *b_t1pfmetPhi;
  TBranch *b_chiso1; 
  TBranch *b_chiso2;  
  TBranch *b_phoiso1; 
  TBranch *b_phoiso2;
  TBranch *b_neuiso1; 
  TBranch *b_neuiso2; 
  TBranch *b_phi1;	  
  TBranch *b_phi2;	  
  TBranch *b_hoe1;	  
  TBranch *b_hoe2;	  
  TBranch *b_sieie1;  
  TBranch *b_sieie2;
  TBranch *b_sceta1;
  TBranch *b_sceta2; 
  TBranch *b_eleveto1;
  TBranch *b_eleveto2;

  // set branch addresses and branch pointers
  treeOrig->SetBranchAddress("run",		&run,		&b_run);
  treeOrig->SetBranchAddress("event",		&event,		&b_event);
  treeOrig->SetBranchAddress("weight",		&weight,	&b_weight);
  treeOrig->SetBranchAddress("nvtx",		&nvtx,		&b_nvtx);
  treeOrig->SetBranchAddress("rho",		&rho,		&b_rho);
  treeOrig->SetBranchAddress("mgg", 		&mgg,		&b_mgg);
  treeOrig->SetBranchAddress("eta1",		&eta1,		&b_eta1); 
  treeOrig->SetBranchAddress("eta2",		&eta2,		&b_eta2); 
  treeOrig->SetBranchAddress("r91",		&r91,		&b_r91);
  treeOrig->SetBranchAddress("r92",		&r92,		&b_r92);
  treeOrig->SetBranchAddress("pt1",		&pt1,		&b_pt1);
  treeOrig->SetBranchAddress("pt2",		&pt2,		&b_pt2);
  treeOrig->SetBranchAddress("ptgg",		&ptgg,		&b_ptgg);
  treeOrig->SetBranchAddress("t1pfmet",		&t1pfmet,	&b_t1pfmet);
  treeOrig->SetBranchAddress("t1pfmetPhi",	&t1pfmetPhi,	&b_t1pfmetPhi);
  treeOrig->SetBranchAddress("chiso1",  	&chiso1,  	&b_chiso1); 
  treeOrig->SetBranchAddress("chiso2",  	&chiso2,  	&b_chiso2);  
  treeOrig->SetBranchAddress("phoiso1", 	&phoiso1, 	&b_phoiso1); 
  treeOrig->SetBranchAddress("phoiso2", 	&phoiso2, 	&b_phoiso2);
  treeOrig->SetBranchAddress("neuiso1", 	&neuiso1, 	&b_neuiso1); 
  treeOrig->SetBranchAddress("neuiso2", 	&neuiso2, 	&b_neuiso2); 
  treeOrig->SetBranchAddress("phi1",    	&phi1,    	&b_phi1);    
  treeOrig->SetBranchAddress("phi2",    	&phi2,    	&b_phi2);    
  treeOrig->SetBranchAddress("hoe1",    	&hoe1,    	&b_hoe1);    
  treeOrig->SetBranchAddress("hoe2",    	&hoe2,    	&b_hoe2);    
  treeOrig->SetBranchAddress("sieie1",  	&sieie1,  	&b_sieie1);  
  treeOrig->SetBranchAddress("sieie2",  	&sieie2,  	&b_sieie2); 
  treeOrig->SetBranchAddress("sceta1",		&sceta1,	&b_sceta1);
  treeOrig->SetBranchAddress("sceta2",		&sceta2,	&b_sceta2);
  treeOrig->SetBranchAddress("eleveto1",	&eleveto1,	&b_eleveto1);
  treeOrig->SetBranchAddress("eleveto2",	&eleveto2,	&b_eleveto2);

  // new variables (needed if variable has diff name in new tree) 
  float mass;
  float leadPt, subleadPt;
  float leadEta, subleadEta;
  float leadScEta, subleadScEta;
  float leadPhoIso, subleadPhoIso;
  float leadChIso, subleadChIso;
  float leadNeutIso, subleadNeutIso; 
  float leadPhi, subleadPhi;
  float leadHoE, subleadHoE;
  float leadSigmaIeIe, subleadSigmaIeIe;
  int   leadPassEleVeto, subleadPassEleVeto;
 
  // setup new trees
  vector<vector<TTree*> > theTreeNew;
  theTreeNew.resize(numMetCat);
  for (UInt_t i=0; i<numMetCat;i++) {
    theTreeNew[i].resize(numPhoCat);
    for (UInt_t j=0; j<numPhoCat;j++) {
      theTreeNew[i][j] = trees[i][j];
      theTreeNew[i][j]->Branch("run",			&run,			"run/I");
      theTreeNew[i][j]->Branch("event",			&event,			"event/I");
      theTreeNew[i][j]->Branch("weight", 		&weight, 		"weight/F");
      theTreeNew[i][j]->Branch("rho",			&rho,			"rho/F");
      theTreeNew[i][j]->Branch("mass", 			&mass, 			"mass/F");
      theTreeNew[i][j]->Branch("nvtx",			&nvtx,			"nvtx/I");
      theTreeNew[i][j]->Branch("leadEta", 		&leadEta, 		"leadEta/F");
      theTreeNew[i][j]->Branch("subleadEta", 		&subleadEta, 		"subleadEta/F");
      theTreeNew[i][j]->Branch("leadPt",		&leadPt,		"leadPt/F");
      theTreeNew[i][j]->Branch("subleadPt",		&subleadPt,		"subleadPt/F");
      theTreeNew[i][j]->Branch("ptgg",			&ptgg,			"ptgg/F");
      theTreeNew[i][j]->Branch("t1pfmet",		&t1pfmet,		"t1pfmet/F");     
      theTreeNew[i][j]->Branch("t1pfmetPhi",		&t1pfmetPhi,		"t1pfmetPhi/F");
      theTreeNew[i][j]->Branch("leadChIso",  		&leadChIso,  		"leadChIso/F"); 
      theTreeNew[i][j]->Branch("subleadChIso",  	&subleadChIso, 		"subleadChIso/F");  
      theTreeNew[i][j]->Branch("leadPhoIso", 		&leadPhoIso, 		"leadPhoIso/F");  
      theTreeNew[i][j]->Branch("subleadPhoIso", 	&subleadPhoIso, 	"subleadPhoIso/F");
      theTreeNew[i][j]->Branch("leadNeutIso", 		&leadNeutIso,	 	"leadNeutIso/F"); 
      theTreeNew[i][j]->Branch("subleadNeutIso",	&subleadNeutIso, 	"subleadNeutIso/F"); 
      theTreeNew[i][j]->Branch("leadPhi",    		&leadPhi,    		"leadPhi/F");    
      theTreeNew[i][j]->Branch("subleadPhi",   		&subleadPhi,    	"subleadPhi/F");    
      theTreeNew[i][j]->Branch("leadHoE",    		&leadHoE,    		"leadHoE/F");    
      theTreeNew[i][j]->Branch("subleadHoE",   		&subleadHoE,    	"subleadHoE/F");    
      theTreeNew[i][j]->Branch("leadSigmaIeIe",		&leadSigmaIeIe,  	"leadSigmaIeIe/F");  
      theTreeNew[i][j]->Branch("subleadSigmaIeIe",  	&subleadSigmaIeIe,  	"subleadSigmaIeIe/F");
      theTreeNew[i][j]->Branch("leadScEta",		&leadScEta,		"leadScEta/F"); 
      theTreeNew[i][j]->Branch("subleadScEta",		&subleadScEta,		"subleadScEta/F"); 
      theTreeNew[i][j]->Branch("leadPassEleVeto",	&leadPassEleVeto,	"leadPassEleVeto/I");
      theTreeNew[i][j]->Branch("subleadPassEleVeto",	&subleadPassEleVeto,	"subleadPassEleVeto/I");
    }// end loop oever new trees in pho cat 
  }// end loop over new trees in met cat

  // make vectors to store if passing different met cuts
  vector<bool> passMet;
  passMet.resize(numMetCat);

  bool EB1, EB2;
  bool EE1, EE2;

  bool inEB, inEE;
  bool hiR9, loR9;

  for (UInt_t i=0; i<nentriesOrig; i++){
    treeOrig->GetEntry(i);

    if (mgg >= 100 && mgg <= 200){

      // split events by eta
      EB1 = false;
      EB2 = false;
      EE1 = false;
      EE2 = false;
      if (fabs(eta1)>=1.5) EE1=true;
      if (fabs(eta2)>=1.5) EE2=true;
      if (fabs(eta1)<1.5)  EB1=true;
      if (fabs(eta2)<1.5)  EB2=true; 
      inEE=false;
      inEB=false;
      if (EB1 && EB2) inEB=true;
      else if (EE1 || EE2) inEE=true;
      
      // split events by r9
      hiR9 = false;
      loR9 = false;
      if (r91 > 0.94 && r92 > 0.94) hiR9 = true;
      else if (r91 > 0.94 || r92 > 0.94) loR9 = true;

      // split events in categories of met
      for (UInt_t met=0; met<numMetCat; met++){
	passMet[met]=false;
        if (met == numMetCat-2 && t1pfmet >= met*metSpacing) passMet[met]=true;			// for 2nd to last bin take all events above met cut
        if (met == numMetCat-1) passMet[met] = true; 						// for last bin take all events NO met cut
        else if (t1pfmet >= met*metSpacing && t1pfmet < (met+1)*metSpacing ) passMet[met]=true; // look at bins in met with value of metSpacing
      }// end met loop

      // set the new variables (i.e. renamed from old tree)
      mass = mgg;
      leadPt = pt1;			subleadPt = pt2;
      leadEta = eta1;			subleadEta = eta2;  
      leadPhi = phi1;			subleadPhi = phi2;
      leadPhoIso = phoiso1;		subleadPhoIso = phoiso2;
      leadNeutIso = neuiso1;		subleadNeutIso = neuiso2;
      leadChIso = chiso1;		subleadChIso = chiso2;
      leadSigmaIeIe = sieie1;		subleadSigmaIeIe = sieie2;
      leadHoE = hoe1;			subleadHoE = hoe2;
      leadPassEleVeto = eleveto1;	subleadPassEleVeto = eleveto2;
      leadScEta = sceta1;		subleadScEta = sceta2;
      
      // fill the trees for events in the different categories
      for (UInt_t met=0; met<numMetCat; met++){
        newDir[met]->cd();
        if (passMet[met]){
            theTreeNew[met][0]->Fill();			 // pho=0 noPhoSel 
	    if (inEB && hiR9) theTreeNew[met][1]->Fill();// pho=1 EBHighR9
	    if (inEB && loR9) theTreeNew[met][2]->Fill();// pho=2 EBLowR9
	    if (inEE && hiR9) theTreeNew[met][3]->Fill();// pho=3 EEHighR9
	    if (inEE && loR9) theTreeNew[met][4]->Fill();// pho=4 EELowR9
        }
      }//end filling loop

    }// end mgg presel 
  }// end loop over entries in orig tree


  // write new file
  fileNew->ls();
  fileNew->cd();
  for (UInt_t met=0; met<numMetCat; met++){
    newDir[met]->cd();
    for (UInt_t pho=0; pho<numPhoCat; pho++){
      theTreeNew[met][pho]->Write();
    }
  }
  
  // close files
  fileNew->Close();
  fileOrig->cd();
  fileOrig->Close();

}// end fitterFormatting
