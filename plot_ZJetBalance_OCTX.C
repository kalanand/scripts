#include "../src/JetUtilMC.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>



TString ZmassWindow = "abs(mZee-91.2)<5.0";




// Z pT binning
const int nZPtBins=8;
const double theZPt[nZPtBins+1] = {25,30,36,43,51,73,104,200,400};


////////////////// Main function: loop over all jet algos ////////

void plot_ZJetBalance_OCTX() {

   for(int i=0; i<6; i++) {plot_ZJetBalance(i);}

}






////////////////// pT balance for a given jet algo ////////



void plot_ZJetBalance( int JetAlgo ) {


  TString etaCut = Form(" && abs(JetRecoEta[%d][0])<1.3", JetAlgo);
  TString dPhiCut =Form(" && JetRecoDphi[%d][0] > 2.94",  JetAlgo);
  TString secondJetCut = Form(" && JetRecoPt[%d][1]/Z_Pt < 0.1", JetAlgo);

  std::string algoSt;
  if(JetAlgo==0) algoSt = "_ic5";
  if(JetAlgo==1) algoSt = "_sc5";
  if(JetAlgo==2) algoSt = "_sc7";
  if(JetAlgo==3) algoSt = "_kt4";
  if(JetAlgo==4) algoSt = "_kt6";
  if(JetAlgo==5) algoSt = "_ak5";

  std::string textFilename = "ZjetPtBalance" + algoSt + ".txt";
  std::string JetResponseFilename 
    = "Histograms_ZjetResponse" + algoSt + ".root";



  ///// open all the input files
  
  TFile* file = TFile::Open("SD_Ele15_Zee.root");
  TTree* tree = (TTree*)file->Get("ZJet");
  TFile* file2 = TFile::Open("SD_Ele15_QCD.root");
  TTree* tree2 = (TTree*)file2->Get("ZJet");
  TFile* file3 = TFile::Open("SD_Ele15_Wenu.root");
  TTree* tree3 = (TTree*)file3->Get("ZJet");
  TFile* file4 = TFile::Open("SD_Ele15_TTbar.root");
  TTree* tree4 = (TTree*)file4->Get("ZJet");



  ////// open the output ROOT file
  TFile respHistFile(JetResponseFilename.c_str(),"RECREATE");
  TString drawSt   = "";
  TString drawSt2  = "";
  TString rGen     = "";
  TString rReco    = "";
  TString pGen     = "";
  TString pReco    = "";
  TString zReco    = "";
  TString ptBin    = "";
  TString zptcut   = "";
  TString allcuts  = "";
  TH1F* responseHistGen(0);
  TH1F* responseHistReco(0);  
  TH1F* genJetpt(0);
  TH1F* caloJetpt(0);
  TH1F* recoZpt(0);


  ////////////////// pT balance for a given Z pT bin ////////

  for(int i=0; i<nZPtBins; i++) {

    ptBin = Form("%d_%d", (int) theZPt[i], (int) theZPt[i+1] );
    rGen  =  TString("responseHistGen_")+ptBin;
    rReco =  TString("responseHistReco_")+ptBin;
    pGen  =  TString("genJetpt_")+ptBin;
    pReco =  TString("caloJetpt_")+ptBin;
    zReco =  TString("recoZpt_")+ptBin;
    zptcut = Form(" && Z_Pt>%f && Z_Pt<%f", theZPt[i], theZPt[i+1] );
    allcuts = ZmassWindow + etaCut + secondJetCut + zptcut;

    responseHistGen  = new TH1F( rGen,"",   10, 0.0, 2.0);
    responseHistReco = new TH1F( rReco, "", 10, 0.0, 2.0);  
    genJetpt         = new TH1F( pGen,"",   20, 0,   400);
    caloJetpt        = new TH1F( pReco, "", 20, 0,   400);
    recoZpt          = new TH1F( zReco,"",  20, 0,   400);
    responseHistReco->Sumw2();
    responseHistGen->Sumw2();
    genJetpt->Sumw2();
    caloJetpt->Sumw2();
    recoZpt->Sumw2();


    drawSt = "Z_Pt>>";
    drawSt2 = "Z_Pt>>+";
    tree->Draw(  drawSt+zReco,  allcuts, "goff");
    tree2->Draw( drawSt2+zReco, allcuts, "goff");
    tree3->Draw( drawSt2+zReco, allcuts, "goff");
    tree4->Draw( drawSt2+zReco, allcuts, "goff");

    drawSt = Form("JetRecoPt[%d][0]>>", JetAlgo);
    drawSt2 = Form("JetRecoPt[%d][0]>>+", JetAlgo);
    tree->Draw(  drawSt+pReco,   allcuts, "goff");
    tree2->Draw( drawSt2+pReco,  allcuts, "goff");
    tree3->Draw( drawSt2+pReco,  allcuts, "goff");
    tree4->Draw( drawSt2+pReco,  allcuts, "goff");


    drawSt = Form("JetRecoPt[%d][0]/Z_Pt>>", JetAlgo);
    drawSt2 = Form("JetRecoPt[%d][0]/Z_Pt>>+", JetAlgo);
    tree->Draw(  drawSt+rReco,  allcuts, "goff");
    tree2->Draw( drawSt2+rReco, allcuts, "goff");
    tree3->Draw( drawSt2+rReco, allcuts, "goff");
    tree4->Draw( drawSt2+rReco, allcuts, "goff");


//     drawSt = Form("JetGenPt[%d][0]>>", JetAlgo);
//     drawSt2 = Form("JetGenPt[%d][0]>>+", JetAlgo);
//     tree->Draw(  drawSt+pGen,   allcuts, "goff");
//     tree2->Draw( drawSt2+pGen, allcuts, "goff");
//     tree3->Draw( drawSt2+pGen, allcuts, "goff");
//     tree4->Draw( drawSt2+pGen, allcuts, "goff");


//     drawSt = Form("JetGenPt[%d][0]/Z_Pt>>", JetAlgo);
//     drawSt2 = Form("JetGenPt[%d][0]/Z_Pt>>+", JetAlgo);
//     tree->Draw(  drawSt+rGen,  allcuts, "goff");
//     tree2->Draw( drawSt2+rGen, allcuts, "goff");
//     tree3->Draw( drawSt2+rGen, allcuts, "goff");
//     tree4->Draw( drawSt2+rGen, allcuts, "goff");


    // Now write all the histograms in a ROOT file
    respHistFile.cd();
    responseHistGen->Write();
    responseHistReco->Write();
    genJetpt->Write();
    caloJetpt->Write();
    recoZpt->Write();


    // clean up the memory
    delete responseHistGen;
    delete responseHistReco;
    delete genJetpt;
    delete caloJetpt;
    delete recoZpt;
  }

  respHistFile.Close();
}





