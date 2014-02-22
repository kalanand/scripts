#include <string>
#include <iostream>
#include <algorithm>
#include <vector>


const char* basedir = "/uscms_data/d2/kalanand/trash/";
//const char* basedir = "./";
const double lumi = 100.0;
// PYTHIA pT-hat binning
const int nGenPtBins=10;
char* pthatBin[nGenPtBins] = 
  {"0_15","15_20","20_30","30_50","50_80","80_120","120_170",
   "170_230","230_300","300_Inf"};

// const char* cmEnergy = "10TeV";
// int nEvents[nGenPtBins] = 
//   { 219730, 220000, 162570, 158606, 110000, 
//     110587, 116496, 131000, 108542, 106355};

// double crosssection[nGenPtBins] = 
//   {  6430.0, 230.0, 211.0, 142.0, 56.8, 18.8, 5.4, 1.55, 0.45, 0.20 };


const char* cmEnergy = "7TeV";

int nEvents[nGenPtBins] = 
  { 214330, 209740, 175590, 191265, 118383, 
    143698, 128045, 138960, 110720, 104900};

double crosssection[nGenPtBins] = 
  { 4434.0, 145.4, 131.8, 84.38, 32.35, 9.981, 2.760, 0.7241, 0.1946, 0.07627 };



// electron id and isolation
const char* eid = "iseMinusLoose && isePlusLoose";
const char* tkIsoCuts = "&& ePlus_trackiso<10.0 && eMinus_trackiso<10.0";
const char* ecalIsoCuts = "&& ePlus_ecaliso<10.0 && eMinus_ecaliso<10.0"; 
const char* hcalIsoCuts = "&& ePlus_hcaliso<6.0 && eMinus_hcaliso<6.0"; 
const char* Zsel = " && abs(mZee-91.2)<10.0";

const char* jetsel = "&& abs(JetCorEta[5][0])<3.0 && JetCorPt[5][0]>20.0 && JetCorEmEnergyFraction[5][0]>0.05 && JetCorEmEnergyFraction[5][0]<0.9"; 
const char* pfsel = "&& abs(JetPFEta[5][0])<3.0 && JetPFPt[5][0]>20.0 && JetPFNeutralEmEnergyFraction[5][0]>0.05 && JetPFNeutralEmEnergyFraction[5][0]<0.9"; 
const std::string genCut = "(abs(JetGenEta[5][0])<3.0 && JetGenPt[5][0]>20.0)";

const std::string electronCuts = std::string(eid)+tkIsoCuts+ecalIsoCuts+hcalIsoCuts;
const std::string allCuts = electronCuts + Zsel+ jetsel;
const std::string allCutsPF  = electronCuts + Zsel+ pfsel;


const char* jetsel2 = "&& abs(JetCorEta[5][1])<3.0 && JetCorPt[5][1]>20.0 && JetCorEmEnergyFraction[5][1]>0.05 && JetCorEmEnergyFraction[5][1]<0.9"; 
const char* pfsel2 = "&& abs(JetPFEta[5][1])<3.0 && JetPFPt[5][1]>20.0 && JetPFNeutralEmEnergyFraction[5][1]>0.05 && JetPFNeutralEmEnergyFraction[5][1]<0.9"; 
const std::string genCut2 = "(abs(JetGenEta[5][1])<3.0 && JetGenPt[5][1]>20.0)";
const std::string allCuts2 = electronCuts + Zsel+ jetsel2;
const std::string allCutsPF2  = electronCuts + Zsel+ pfsel2;


const char* jetsel3 = "&& abs(JetCorEta[5][2])<3.0 && JetCorPt[5][2]>20.0 && JetCorEmEnergyFraction[5][2]>0.05 && JetCorEmEnergyFraction[5][2]<0.9"; 
const char* pfsel3 = "&& abs(JetPFEta[5][2])<3.0 && JetPFPt[5][2]>20.0 && JetPFNeutralEmEnergyFraction[5][2]>0.05 && JetPFNeutralEmEnergyFraction[5][2]<0.9"; 
const std::string genCut3 = "(abs(JetGenEta[5][2])<3.0 && JetGenPt[5][2]>20.0)";
const std::string allCuts3 = electronCuts + Zsel+ jetsel3;
const std::string allCutsPF3  = electronCuts + Zsel+ pfsel3;




// some default hotogram binnings and title
const double MIN = 0.0;
const double MAX = 400.0;
const int BINs = 20;


// declare all histograms
TH1D* ZPt_Gen;
TH1D* ZPt_Reco;
TH1D* GenJetPt;
TH1D* CaloJetPt;
TH1D* PFJetPt;
TH1D* emPt_Gen;
TH1D* emPt_Reco;
TH1D* epPt_Gen;
TH1D* epPt_Reco;
TH1D* Zmass;

TH1D* ZEta_Gen;
TH1D* ZEta_Reco;
TH1D* GenJetEta;
TH1D* CaloJetEta;
TH1D* PFJetEta;

TH1D* GenJetPt2;
TH1D* CaloJetPt2;
TH1D* PFJetPt2;

TH1D* GenJetPt3;
TH1D* CaloJetPt3;
TH1D* PFJetPt3;

TH1D* CaloJetPtRatio2over1;
TH1D* CaloJetPtRatio3over2;
TH1D* PFJetPtRatio2over1;
TH1D* PFJetPtRatio3over2;
TH1D* GenJetPtRatio2over1;
TH1D* GenJetPtRatio3over2;


TH1D* nJetsCalo;
TH1D* nJetsPF;
TH1D* nJetsGen;

/////////////////////////////////////////////////////

void plotSpectrum() {
  CreateHistograms();
  char* append = "";

  // ****** get the weight value for each pthat bin ******* //
  for(int i=0; i<nGenPtBins; i++) {
    // for(int i=5; i<6; i++) {
    if(i!=0) append = "+";
    DrawTreeForAllHistos( i, append );
  }


  ///////////////////////////////////////

  nJetsCalo->SetBinContent( 1, Zmass->Integral() ); 
  nJetsCalo->SetBinError( 1, sqrt(Zmass->Integral()) );
  nJetsPF->SetBinContent( 1, Zmass->Integral());
  nJetsPF->SetBinError( 1, sqrt(Zmass->Integral()) );
  nJetsGen->SetBinContent( 1, Zmass->Integral());
  nJetsGen->SetBinError( 1, sqrt(Zmass->Integral()) );

  nJetsCalo->SetBinContent( 2, CaloJetPt->Integral() ); 
  nJetsCalo->SetBinError( 2, sqrt(CaloJetPt->Integral()) );
  nJetsPF->SetBinContent( 2, PFJetPt->Integral());
  nJetsPF->SetBinError( 2, sqrt(PFJetPt->Integral()) );
  nJetsGen->SetBinContent( 2, GenJetPt->Integral());
  nJetsGen->SetBinError( 2, sqrt(GenJetPt->Integral()) );

  nJetsCalo->SetBinContent( 3, CaloJetPt2->Integral()); 
  nJetsCalo->SetBinError( 3, sqrt(CaloJetPt2->Integral()) );
  nJetsPF->SetBinContent( 3, PFJetPt2->Integral());
  nJetsPF->SetBinError( 3, sqrt(PFJetPt2->Integral()) );
  nJetsGen->SetBinContent( 3, GenJetPt2->Integral());
  nJetsGen->SetBinError( 3, sqrt(GenJetPt2->Integral()) );

  nJetsCalo->SetBinContent( 4, CaloJetPt3->Integral()); 
  nJetsCalo->SetBinError( 4, sqrt(CaloJetPt3->Integral()) );
  nJetsPF->SetBinContent( 4, PFJetPt3->Integral());
  nJetsPF->SetBinError( 4, sqrt(PFJetPt3->Integral()) );
  nJetsGen->SetBinContent( 4, GenJetPt3->Integral());
  nJetsGen->SetBinError( 4, sqrt(GenJetPt3->Integral()) );


  /////////////////////////////////////
  ZPt_Gen->Sumw2();
  ZPt_Reco->Sumw2();
  GenJetPt->Sumw2();
  CaloJetPt->Sumw2();
  PFJetPt->Sumw2();

  emPt_Gen->Sumw2();
  emPt_Reco->Sumw2();
  epPt_Gen->Sumw2();
  epPt_Reco->Sumw2();
  Zmass->Sumw2();
  ZEta_Gen->Sumw2();
  ZEta_Reco->Sumw2();
  GenJetEta->Sumw2();
  CaloJetEta->Sumw2();
  PFJetEta->Sumw2();

  GenJetPt2->Sumw2();
  CaloJetPt2->Sumw2();
  PFJetPt2->Sumw2();

  GenJetPt3->Sumw2();
  CaloJetPt3->Sumw2();
  PFJetPt3->Sumw2();


  // set proper style for plots
  gROOT->ProcessLine(".L mystyle.C");
  setTDRStyle();
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPadLeftMargin(0.18);
  tdrStyle->SetPadRightMargin(0.10);
  tdrStyle->SetPadBottomMargin(0.16);
  tdrStyle->SetLegendBorderSize(0);
  tdrStyle->SetTitleYOffset(1.5);

  makeplotTwo(*ZPt_Reco, *ZPt_Gen, "ZPt_spectrum", 1);
  // makeplotTwo( *emPt_Reco, *emPt_Gen, "emPt_spectrum", 0);
  // makeplotTwo( *epPt_Reco, *epPt_Gen, "epPt_spectrum", 0);
  // makeplotTwo(*ZEta_Reco, *ZEta_Gen, "Z_eta_spectrum", 2);

  makeplotThree(*CaloJetPt, *GenJetPt, *PFJetPt, "Jet_spectrum", 1);
  makeplotThree(*CaloJetPt2, *GenJetPt2, *PFJetPt2, "Jet_spectrum2", 1);
  makeplotThree(*CaloJetPt3, *GenJetPt3, *PFJetPt3, "Jet_spectrum3", 1);

  makeplotThree(*CaloJetEta, *GenJetEta, *PFJetEta, "Jet_eta_spectrum", 2);
  //makeplot( *Zmass, "Zmass", "e", 0);



  //make ratio histograms
  CaloJetPtRatio2over1 = makeRatioHist( *CaloJetPt, *CaloJetPt2, "CaloJetPtRatio2over1",150);
  CaloJetPtRatio3over2 = makeRatioHist( *CaloJetPt2, *CaloJetPt3, "CaloJetPtRatio2over1",100);
  PFJetPtRatio2over1 = makeRatioHist( *PFJetPt, *PFJetPt2, "PFJetPtRatio2over1",150);
  PFJetPtRatio3over2 = makeRatioHist( *PFJetPt2, *PFJetPt3, "PFJetPtRatio3over2",100);
  GenJetPtRatio2over1 = makeRatioHist( *GenJetPt, *GenJetPt2, "GenJetPtRatio2over1",150);
  GenJetPtRatio3over2 = makeRatioHist( *GenJetPt2, *GenJetPt3, "GenJetPtRatio3over2",100);
  CaloJetPtRatio2over1->GetYaxis()->SetTitle("n^{Jet}_{#geq 2} / n^{Jet}_{#geq 1}");
  CaloJetPtRatio3over2->GetYaxis()->SetTitle("n^{Jet}_{#geq 3} / n^{Jet}_{#geq 2}");

  makeplotThree(*CaloJetPtRatio2over1, *GenJetPtRatio2over1, *PFJetPtRatio2over1, "JetPtRatio2over1", 0);
  makeplotThree(*CaloJetPtRatio3over2, *GenJetPtRatio3over2, *PFJetPtRatio3over2, "JetPtRatio3over2", 0);
  makeplotThree( *nJetsCalo, *nJetsGen, *nJetsPF, "JetMultiplicity", 1);

  //////////////////////////////////////
  //tdrStyle->SetPadLeftMargin(0.22);
  TFile outfile( (std::string("histograms_")+cmEnergy+".root").c_str(), "RECREATE"  );
  outfile.cd();
  ZPt_Gen->Write();
  ZPt_Reco->Write();
  GenJetPt->Write();
  CaloJetPt->Write();
  PFJetPt->Write();
  emPt_Gen->Write();
  emPt_Reco->Write();
  epPt_Gen->Write();
  epPt_Reco->Write();
  Zmass->Write();

  ZEta_Gen->Write();
  ZEta_Reco->Write();
  GenJetEta->Write();
  CaloJetEta->Write();
  PFJetEta->Write();

  GenJetPt2->Write();
  CaloJetPt2->Write();
  PFJetPt2->Write();

  GenJetPt3->Write();
  CaloJetPt3->Write();
  PFJetPt3->Write();

  nJetsCalo->Write();
  nJetsPF->Write();
  nJetsGen->Write();

  CaloJetPtRatio2over1->Write();
  CaloJetPtRatio3over2->Write();
  PFJetPtRatio2over1->Write();
  PFJetPtRatio3over2->Write();
  GenJetPtRatio2over1->Write();
  GenJetPtRatio3over2->Write();

  outfile.Close();
}




void DrawTreeForAllHistos( int fileindex, char* append) {

  std::string filename = std::string(basedir);
  filename.append("Summer09-");
  filename.append(cmEnergy);
  filename.append("-ZeeJets-Pt_");
  filename.append( pthatBin[fileindex] );
  filename.append( ".root" );
  
  TFile f( filename.c_str(), "read");
  TTree* atree = (TTree*) f.Get("ZJet");
  double weight = lumi * crosssection[fileindex] / nEvents[fileindex];
  atree->SetWeight( weight );

  DrawTree( *atree, *Zmass,      "mZee>>", electronCuts, append );
  DrawTree( *atree, *emPt_Reco,  "eMinusPt>>", allCuts,         append );
  DrawTree( *atree, *epPt_Reco,  "ePlusPt>>",  allCuts,         append );
  DrawTree( *atree, *ZPt_Reco,   "Z_Pt>>",  allCuts,            append );
  DrawTree( *atree, *ZEta_Reco,  "Z_Eta>>",   allCuts,          append );
  DrawTree( *atree, *CaloJetPt,  "JetCorPt[5][0]>>",  allCuts,  append );
  DrawTree( *atree, *CaloJetEta, "JetCorEta[5][0]>>", allCuts,  append );
  DrawTree( *atree, *PFJetPt,  "JetPFPt[5][0]>>",  allCutsPF,  append );
  DrawTree( *atree, *PFJetEta, "JetPFEta[5][0]>>", allCutsPF,  append );   
 
  // Now plot the generator level quantities
  DrawTree( *atree, *ZPt_Gen,    "Z_PtGen>>",         genCut,  append );
  DrawTree( *atree, *ZEta_Gen,   "Z_EtaGen>>",        genCut,  append );
  DrawTree( *atree, *emPt_Gen,   "eMinusPtGen>>",     genCut,  append );
  DrawTree( *atree, *epPt_Gen,   "ePlusPtGen>>",      genCut,  append );
  DrawTree( *atree, *GenJetPt,   "JetGenPt[5][0]>>",  genCut,  append );   
  DrawTree( *atree, *GenJetEta,  "JetGenEta[5][0]>>", genCut,  append );

  DrawTree( *atree, *CaloJetPt2,  "JetCorPt[5][1]>>",  allCuts2,  append );
  DrawTree( *atree, *PFJetPt2,  "JetPFPt[5][1]>>",  allCutsPF2,  append );
  DrawTree( *atree, *GenJetPt2,   "JetGenPt[5][1]>>",  genCut2,  append );   

  DrawTree( *atree, *CaloJetPt3,  "JetCorPt[5][2]>>",  allCuts3,  append );
  DrawTree( *atree, *PFJetPt3,  "JetPFPt[5][2]>>",  allCutsPF3,  append );
  DrawTree( *atree, *GenJetPt3,   "JetGenPt[5][2]>>",  genCut3,  append );   
} 




void DrawTree( TTree& tree, TH1& hist, const char* branchname, 
               const std::string condition, const char* append = "") {

  gROOT->cd();

  std::string plot = std::string(branchname);
  plot.append(append);
  plot.append(hist.GetName());

  tree.Draw( plot.c_str(), condition.c_str(), "goff" );
}





void makeplotThree(TH1& hist1, TH1& hist2, TH1& hist3, const char* plotname, int log) {

  hist1.SetLineColor(4);
  hist1.SetMarkerColor(4);
  hist3.SetLineColor(2);
  hist3.SetMarkerColor(2);

  gStyle->SetOptStat(0);
  TCanvas* can = new TCanvas( plotname, plotname, 500, 500);
  hist1.Draw( );
  hist2.Draw( "same" );
  hist3.Draw( "same" );
  TLegend *leg = new TLegend(0.48,0.7,0.89,0.92);
  leg->AddEntry( &hist1,"CaloJet","LP");
  leg->AddEntry( &hist3,"PF Jet","LP");
  leg->AddEntry( &hist2,"GenJet","L");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
  can->SetLogy( log );
  if(log==2) {
    can->SetLogy( 1 );
    hist1.GetYaxis()->SetMoreLogLabels();
  }

  std::string plot(cmEnergy);
  plot.append(plotname);
  can->SaveAs( (plot+".eps").c_str() );
  can->SaveAs( (plot+".gif").c_str() );
  can->SaveAs( (plot+".root").c_str() );
  // delete can;
  cout << hist1.Integral() << endl;
}




void makeplotTwo(TH1& hist1, TH1& hist2, const char* plotname, int logy) {

  hist1.SetLineColor(4);
  hist1.SetMarkerColor(4);
  gStyle->SetOptStat(0);
  TCanvas* can = new TCanvas( plotname, plotname, 500, 500);
  hist1.Draw( );
  hist2.Draw( "same" );
  TLegend *leg = new TLegend(0.55,0.8,0.89,0.92);
  leg->AddEntry( &hist1,"Reconstructed","LP");
  leg->AddEntry( &hist2,"Generated","L");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
  can->SetLogy( logy );
  if(logy==2) hist1.GetYaxis()->SetMoreLogLabels();

  std::string plot(cmEnergy);
  plot.append(plotname);
  can->SaveAs( (plot+".eps").c_str() );
  can->SaveAs( (plot+".gif").c_str() );
  can->SaveAs( (plot+".root").c_str() );
  // delete can;
  cout << hist1.Integral() << endl;
}



void makeplot(TH1& hist, const char* plotname, const char* option, int log) {

  gStyle->SetOptStat(0);
  TCanvas* can = new TCanvas( plotname, plotname, 500, 500);
  hist.Draw( option );
  can->SetLogy( log );

  std::string plot(cmEnergy);
  plot.append(plotname);
  can->SaveAs( (plot+".eps").c_str() );
  can->SaveAs( (plot+".gif").c_str() );
  can->SaveAs( (plot+".root").c_str() ); 
  // delete can;
  cout << hist.Integral() << endl;
}





void SetAxisLabels(TH1& hist, char* xtitle, char* ytitle="",
		   double xoffset=1.1, double yoffset=1.4) {

  TAxis* x = hist.GetXaxis();
  TAxis* y = hist.GetYaxis();
  x->SetTitle(xtitle);
  x->SetTitleSize(0.06);
  x->SetLabelSize(0.05);
  x->SetTitleOffset(xoffset);
  x->SetNdivisions(505);
  y->SetTitle(ytitle);
  y->SetTitleSize(0.06);
  y->SetLabelSize(0.05);
  y->SetTitleOffset(yoffset);
  y->SetNoExponent();
  hist.SetLineWidth(2);
  hist.SetMarkerStyle(20);

  std::stringstream str;
  str << "Events / " << (int) lumi << " pb^{-1}   ";
  std::string  defYtitle = str.str();
  if(ytitle=="") y->SetTitle( defYtitle.c_str() );
}





void CreateHistograms() {

  ZPt_Gen = new TH1D("ZPt_Gen","", BINs, MIN, MAX);
  ZPt_Reco = new TH1D("ZPt_Reco","", BINs, MIN, MAX);
  GenJetPt = new TH1D("GenJetPt","", BINs, MIN, MAX);
  CaloJetPt = new TH1D("CaloJetPt","", BINs, MIN, MAX);
  PFJetPt = new TH1D("PFJetPt","", BINs, MIN, MAX);

  emPt_Gen = new TH1D("emPt_Gen","", 100, 0, 100);
  emPt_Reco = new TH1D("emPt_Reco","", 100, 0, 100);
  epPt_Gen = new TH1D("epPt_Gen","", 100, 0, 100);
  epPt_Reco = new TH1D("epPt_Reco","", 100, 0, 100);
  Zmass = new TH1D("Zmass","", 120, 60, 120);

  ZEta_Gen   = new TH1D("ZEta_Gen","",   30, -3, 3);
  ZEta_Reco  = new TH1D("ZEta_Reco","",  30, -3, 3);
  GenJetEta  = new TH1D("GenJetEta","",  30, -3, 3);
  CaloJetEta = new TH1D("CaloJetEta","", 30, -3, 3);
  PFJetEta = new TH1D("PFJetEta","",     30, -3, 3);

  GenJetPt2 = new TH1D("GenJetPt2","", BINs, MIN, MAX);
  CaloJetPt2 = new TH1D("CaloJetPt2","", BINs, MIN, MAX);
  PFJetPt2 = new TH1D("PFJetPt2","", BINs, MIN, MAX);

  GenJetPt3 = new TH1D("GenJetPt3","", BINs, MIN, MAX);
  CaloJetPt3 = new TH1D("CaloJetPt3","", BINs, MIN, MAX);
  PFJetPt3 = new TH1D("PFJetPt3","", BINs, MIN, MAX);

  nJetsCalo = new TH1D("nJetsCalo", "", 4, 0, 4);
  nJetsPF = new TH1D("nJetsPF", "",     4, 0, 4);
  nJetsGen = new TH1D("nJetsGen", "",   4, 0, 4);

  SetAxisLabels( *ZPt_Gen, "Z p_{T} (GeV/c)");
  SetAxisLabels( *ZPt_Reco, "Z p_{T} (GeV/c)");
  SetAxisLabels( *GenJetPt, "leading jet p_{T} (GeV/c)");
  SetAxisLabels( *CaloJetPt, "leading jet p_{T} (GeV/c)");
  SetAxisLabels( *PFJetPt, "leading jet p_{T} (GeV/c)");


  SetAxisLabels( *emPt_Gen, "e^{-} p_{T} (GeV/c)");
  SetAxisLabels( *emPt_Reco, "e^{-} p_{T} (GeV/c)");
  SetAxisLabels( *epPt_Gen, "e^{+} p_{T} (GeV/c)");
  SetAxisLabels( *epPt_Reco, "e^{+} p_{T} (GeV/c)");
  SetAxisLabels( *Zmass, "m_{e^{+}e^{-}} (GeV/c^{2})");

  SetAxisLabels( *ZEta_Gen,    "#eta");
  SetAxisLabels( *ZEta_Reco,   "#eta");
  SetAxisLabels( *GenJetEta,   "#eta");
  SetAxisLabels( *CaloJetEta,  "#eta");
  SetAxisLabels( *PFJetEta,  "#eta");

  SetAxisLabels( *GenJetPt2, "second jet p_{T} (GeV/c)");
  SetAxisLabels( *CaloJetPt2, "second jet p_{T} (GeV/c)");
  SetAxisLabels( *PFJetPt2, "second jet p_{T} (GeV/c)");

  SetAxisLabels( *GenJetPt3, "third jet p_{T} (GeV/c)");
  SetAxisLabels( *CaloJetPt3, "third jet p_{T} (GeV/c)");
  SetAxisLabels( *PFJetPt3, "third jet p_{T} (GeV/c)");

  SetAxisLabels( *nJetsCalo, "Jet multiplicity, n^{Jets}_{#geq}");
  SetAxisLabels( *nJetsPF, "Jet multiplicity, n^{Jets}_{#geq}");
  SetAxisLabels( *nJetsGen, "Jet multiplicity, n^{Jets}_{#geq}");

  ZEta_Reco->GetXaxis()->SetNdivisions(510);
  CaloJetEta->GetXaxis()->SetNdivisions(510);
 
}



// returns  (h2)/(h1-h2)

TH1D* makeRatioHist(TH1D& Hist1, TH1D& Hist2, char* name="CaloJetPtRatio2over1", double maxX=100.0) {
  
  TH1D* hist = Hist2.Clone(name);
  TH1D* temphist = Hist1.Clone("temphist");
  //temphist->Add( &Hist2, -1);
  hist->Divide(temphist);
  hist->GetXaxis()->SetTitle("jet p_{T} (GeV/c)");
  hist->SetMinimum(0.0);
  if(maxX > 110.0) hist->SetMaximum(0.35);
  else hist->SetMaximum(0.3);
  hist->GetXaxis()->SetRangeUser(0,maxX);
  delete temphist;
  return hist;
}
