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
void sigModelFit(RooWorkspace*);
void drawPlots(RooWorkspace*, TString, Int_t, Float_t, Float_t, TString);
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

  RooDataSet* signal[nMetCat][nPhoCat];

  TFile* inFile = TFile::Open(Form("%s_new.root",name.Data()));
  if (inFile == (TFile*) NULL) std::cout<<" NOT A VALID FILE " << std::endl;
  
  TTree* sigTree1 = new TTree();

  for (UInt_t met=0; met<nMetCat; met++){
    for (UInt_t pho=0; pho<nPhoCat; pho++){

      sigTree1 = (TTree*)inFile->Get(Form("%s/%s",MetCat[met].Data(),PhoCat[pho].Data()));
      if (sigTree1 == (TTree*) NULL) std::cout << " NOT A VALID TREE " << std::endl;
      sigTree1->SetTitle(name);
      sigTree1->SetName(name);

      signal[met][pho] = new RooDataSet("signal","dataset",sigTree1,*ntplVars,mainCut);//,"weight");
      signal[met][pho]->Print("v");
 
      w->import(*signal[met][pho],Rename(TString::Format("%s_%s",MetCat[met].Data(),PhoCat[pho].Data()))); 

    }// end loop over pho cat
  }// end loop over met cat
  //w->Print("v");
}

void sigModelFit(RooWorkspace* w){

  Float_t mass=125.;
  Float_t minMassFit(mass*0.8);
  Float_t maxMassFit(mass*1.2);

  for (UInt_t met=0; met < nMetCat; met++){
    for (UInt_t pho=0; pho < nPhoCat; pho++){
      

    }
  }

}


void drawPlots(RooWorkspace* w, TString variable, int BINS, float MIN, float MAX, TString mass){
  TString inDir = "";
  TString name = "";
  TString datasetName = "";
  std::vector<TString> MetCat = defineMet(nMetCat);
  std::vector<TString> PhoCat = definePho(nPhoCat);
 
  TFile* f = new TFile("signalPlots.root","RECREATE");
  f->cd();   

  TCanvas* c[nPhoCat]; 
  for (UInt_t pho=0; pho<nPhoCat; pho++){
   c[pho] = new TCanvas(Form("c%s",PhoCat[pho].Data()),"c",1);
  }

  Color_t colorMetCat[nMetCat];
  colorMetCat[0]=kGreen-9;
  colorMetCat[1]=kYellow+8;
  colorMetCat[2]=kTeal-1;
  colorMetCat[3]=kGreen;
  colorMetCat[4]=kTeal;
  colorMetCat[5]=kBlue;
  colorMetCat[6]=kBlue+2;
  colorMetCat[7]=kViolet;
  colorMetCat[8]=kMagenta;
  colorMetCat[9]=kPink-2;
  colorMetCat[10]=kPink+6;
  colorMetCat[11]=kBlack;

  RooDataSet* sigDataSet[nMetCat][nPhoCat];
  RooPlot* sigth1f[nPhoCat];

  for (UInt_t met=0; met<nMetCat; met++){
    for (UInt_t pho=0; pho<nPhoCat; pho++){
      datasetName = TString::Format("%s_%s",MetCat[met].Data(),PhoCat[pho].Data());
      name = TString::Format("h_%s_%s",MetCat[met].Data(),PhoCat[pho].Data()); 

      sigth1f[pho] = w->var("mass")->frame(Range(MIN,MAX),Bins(BINS));
      if(sigth1f[pho]== (RooPlot*) NULL) std::cout<<"VARIABLE NOT FOUND" << std::endl;
      sigDataSet[met][pho] = (RooDataSet*) w->data(datasetName);
    }
  } 

  for (UInt_t pho=0; pho<nPhoCat; pho++){
     c[pho]->cd();
     for (UInt_t met=0; met<nMetCat; met++){
       // plot all met bins in same pho plot
       sigDataSet[met][pho]->plotOn(sigth1f[pho],LineColor(colorMetCat[met]),DrawOption("L"),LineStyle(1),MarkerStyle(0),XErrorSize(0),DataError(RooAbsData::None));
     }
     sigth1f[pho]->Draw(); 
     TLegend* leg = new TLegend(0.55,0.6,0.87,0.88,(TString::Format("%s",PhoCat[pho].Data())),"brNDC");
     for (UInt_t met=0; met<nMetCat; met++){
       leg->AddEntry(sigth1f[pho]->getObject(met),MetCat[met].Data(),"L");
     }
     leg->Draw();
     c[pho]->SetLogy(0);
     c[pho]->SaveAs(TString::Format("plots/%s.png",PhoCat[pho].Data()));
     c[pho]->SetLogy(1);
     c[pho]->SaveAs(TString::Format("plots/%s_log.png",PhoCat[pho].Data()));
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
  std::cout << "Starting SigModelFit" << std::endl;
  sigModelFit(w);

  std::cout << "Making Plots" << std::endl;
  drawPlots(w,"mgg",40,110.,150.,"600");

  w->writeToFile("signals_wkspace.root");
 
}






