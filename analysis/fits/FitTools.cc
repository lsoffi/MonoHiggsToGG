using namespace RooFit;
using namespace RooStats;

static const Int_t nMetCat=12;
static const Int_t nPhoCat=5;
Int_t MinMass = 100;
Int_t MaxMass = 200;

std::vector<TString> defineMet(Int_t metCat){
  std::vector<TString> MetCat;
  MetCat.resize(metCat);
  MetCat[0]="met0_100";
  MetCat[1]="met100_200";
  MetCat[2]="met200_300";
  MetCat[3]="met300_400";
  MetCat[4]="met400_500";
  MetCat[5]="met500_600";
  MetCat[6]="met600_700";
  MetCat[7]="met700_800";
  MetCat[8]="met800_900";
  MetCat[9]="met900_1000";
  MetCat[10]="met1000_1100";
  MetCat[11]="met_all";
  
  return MetCat;
}

std::vector<TString> definePho(Int_t phoCat){

  std::vector<TString> PhoCat;
  PhoCat.resize(phoCat);
  PhoCat[0]="sig_2HDM_mZP600_mA0300_all";
  PhoCat[1]="sig_2HDM_mZP600_mA0300_EBHighR9";
  PhoCat[2]="sig_2HDM_mZP600_mA0300_EBLowR9";
  PhoCat[3]="sig_2HDM_mZP600_mA0300_EEHighR9";
  PhoCat[4]="sig_2HDM_mZP600_mA0300_EELowR9";
  
  return PhoCat;
}

void AddSigData(RooWorkspace*, Float_t);
RooArgSet* defineVariables();


RooArgSet* defineVariables(){

  RooRealVar* mass	= new RooRealVar("mass","m(gg)",100,200,"GeV");
  RooRealVar* leadEta	= new RooRealVar("leadEta","eta(g1)",-10,10,"");
  RooRealVar* subleadEta = new RooRealVar("subleadEta","eta(g2)",-10,10,"");
  RooRealVar* leadR9	= new RooRealVar("leadR9","r9(g1)",-10,10,"");
  RooRealVar* subleadR9	= new RooRealVar("subleadR9","r9(g2)",-10,10,"");
  RooRealVar* nvtx    = new RooRealVar("nvtx","nvtx",0,60,"");
  RooRealVar* weight  = new RooRealVar("weight","weight",-5,5,"");
  RooRealVar* t1pfmet = new RooRealVar("t1pfmet","t1pfmet",0,1200,""); 
  RooRealVar* passHlt = new RooRealVar("passHlt","passHlt",-0.5,1.5,"");

  RooArgSet* ntplVars = new RooArgSet(*mass,*leadEta,*subleadEta,*leadR9,*subleadR9,*nvtx,*t1pfmet,*weight,*passHlt);

  return ntplVars;
}

void AddSigData(RooWorkspace* w, TString Mass){
  TString name = TString::Format("2HDM_mZP%s",Mass.Data());
 
  // Variables
  RooArgSet* ntplVars = defineVariables();
  std::vector<TString> MetCat = defineMet(nMetCat);
  std::vector<TString> PhoCat = definePho(nPhoCat);

  TString mainCut = "mass>= 100 && mass <= 200 && passHlt==1";

  TChain* signalTree[nMetCat][nPhoCat];
  RooDataSet* signal[nMetCat][nPhoCat];

  for (UInt_t met=0; met<nMetCat; met++){
    for (UInt_t pho=0; pho<nPhoCat; pho++){
      signalTree[met][pho] = new TChain();
      signalTree[met][pho]->Add(Form("%s_new.root/%s/%s",name.Data(),MetCat[met].Data(),PhoCat[pho].Data()));
 

      TTree* sigTree1 = new TTree();
      sigTree1->Clone(Form("%s_new.root/%s/%s",name.Data(),MetCat[met].Data(),PhoCat[pho].Data()));
      sigTree1->SetTitle(name);
      sigTree1->SetName(name);


      // RooDataSet from TChain does NOT work
      // Also workspace->import does NOT work for RooRealVar or RooDataSet


      signal[met][pho] = new RooDataSet("signal","dataset",sigTree1,*ntplVars,mainCut,"weight");
      //RooDataSet* signal("signal","dataset",signalTree[met][pho],*ntplVars,mainCut);
      signal[met][pho]->Print("v");
 
      w->import(*signal[met][pho],Rename(TString::Format("Sig_%s_%s_%s",Mass.Data(),MetCat[met].Data(),PhoCat[pho].Data()))); 
      //w->Print("v");

    }
  }



}

void runfits(){

  TString card_name("MonoHiggs.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();

  Float_t Lum = 1263.0; 
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi);
  //w->Print();
 
  std::cout << "Adding Signal Samples" << std::endl;
  AddSigData(w,"600");
  //AddSigData(w,"800");
  //AddSigData(w,"1000");
  //AddSigData(w,"1200");
  //AddSigData(w,"1400");

}






