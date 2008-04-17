#include "../src/JetUtilMC.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>


const bool usingCorrectedCaloJetPt = false;
const bool applyWeights = true;
const double lumi = 100.0;
bool ScaleToLuminosity = true;
TString basedir = "/uscms_data/d1/kalanand/trash/";

// Z pT binning
const double MINPTCUT = 1.0;
const int nZPtBins=15;
const double theZPt[nZPtBins+1] = {25,30,36,43,51,61,73,87,104,124,148,
				   177,212,254,304,364};


Float_t pt[nZPtBins];
Float_t errpt[nZPtBins];
TString ptBin[nZPtBins];

// PYTHIA pT-hat binning
const int nGen=10;
double weight[nGen];
// int nEvents[nGen] = { 215061, 201960, 212760, 179490, 104160, 
// 		      103878, 105112, 105941, 104940, 105376};
int nEvents[nGen];

double crosssection[nGen];
crosssection[0]  =  6430.0;
crosssection[1]  =  230.0 ;
crosssection[2]  =  211.0 ;
crosssection[3]  =  142.0 ;
crosssection[4]  =  56.8  ;
crosssection[5]  =  18.8  ;
crosssection[6]  =  5.4   ;
crosssection[7]  =  1.55  ;
crosssection[8]  =  0.45  ;
crosssection[9]  =  0.20  ;

TString inFileNames[nGen] = {"redigi-Summer08-ZeeJets_Pt_0_15.root",
			     "redigi-Summer08-ZeeJets_Pt_15_20.root",
			     "redigi-Summer08-ZeeJets_Pt_20_30.root",
			     "redigi-Summer08-ZeeJets_Pt_30_50.root",
			     "redigi-Summer08-ZeeJets_Pt_50_80.root",
			     "redigi-Summer08-ZeeJets_Pt_80_120.root",
			     "redigi-Summer08-ZeeJets_Pt_120_170.root",
			     "redigi-Summer08-ZeeJets_Pt_170_230.root",
			     "redigi-Summer08-ZeeJets_Pt_230_300.root",
			     "redigi-Summer08-ZeeJets_Pt_300_Inf.root"};



// TString inFileNames[nGen] = {"Summer08_ZeeJets_Pt_0_15.root",
// 			     "Summer08_ZeeJets_Pt_15_20.root",
// 			     "Summer08_ZeeJets_Pt_20_30.root",
// 			     "Summer08_ZeeJets_Pt_30_50.root",
// 			     "Summer08_ZeeJets_Pt_50_80.root",
// 			     "Summer08_ZeeJets_Pt_80_120.root",
// 			     "Summer08_ZeeJets_Pt_120_170.root",
// 			     "Summer08_ZeeJets_Pt_170_230.root",
// 			     "Summer08_ZeeJets_Pt_230_300.root",
// 			     "Summer08_ZeeJets_Pt_300_Inf.root"};


// constants and user options
const double M_PI = TMath::Pi();
const int NUM_JET_MAX = 10;
const bool storeResponsHistogramsInRootFile = true;
const bool makeplot_ZptBalance = false;
const bool makeplot_Response = false;
const bool storeResponseInTextFile = false;


// set proper style for plots
gROOT->ProcessLine(".L mystyle.C");
setTDRStyle();
tdrStyle->SetErrorX(0.5);
tdrStyle->SetPadLeftMargin(0.14);
tdrStyle->SetPadRightMargin(0.10);
tdrStyle->SetLegendBorderSize(0);
tdrStyle->SetTitleYOffset(1.5);




/////////////////////////////////////////////////////
///////////  Setup the tree branches //////////////

TChain* mychain = new TChain("ZJet"); 
int numFiles = 
//mychain->Add("/uscms_data/d1/kalanand/trash/Summer08_ZeeJets_Pt_*.root");
mychain->Add("/uscms_data/d1/kalanand/trash/redigi-Summer08-ZeeJets_Pt_*.root");
Float_t JetRecoPt[10][10];
Float_t JetCorPt[10][10];
Float_t JetCorEta[10][10];
Float_t JetCorPhi[10][10];
Float_t JetGenPt[10][10];
Float_t JetRecoEta[10][10];
Float_t JetGenEta[10][10];
Float_t JetRecoPhi[10][10];
Float_t JetGenPhi[10][10];
Float_t JetRecoType[10][10];
Float_t Z_Pt;
Float_t Z_Phi;
Float_t Z_Eta;
Float_t Z_PtGen;
Float_t Z_PhiGen;
Float_t ePlusPt;
Float_t eMinusPt;
Float_t ePlusPtGen;
Float_t eMinusPtGen;
Float_t ePlusEta;
Float_t eMinusEta;
Float_t ePlusEtaGen;
Float_t eMinusEtaGen;
Float_t ePlusPhi;
Float_t eMinusPhi;
Float_t ePlusPhiGen;
Float_t eMinusPhiGen;
Float_t mZeeGen;
Float_t mZee;
Float_t ePlusE;
Float_t eMinusE;
Float_t ePlus_ecaliso;
Float_t eMinus_ecaliso;
Float_t ePlus_hcaliso;
Float_t eMinus_hcaliso;
Float_t ePlus_trackiso;	   
Float_t eMinus_trackiso;
Int_t NumRecoJets;
Int_t NumGenJets;  

Bool_t iseMinusLoose;
Bool_t isePlusLoose;
   
float R1[NUM_JET_MAX] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
float R2[NUM_JET_MAX] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
float p[NUM_JET_MAX] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  


////////////////// Main function: loop over all jet algos ////////

void plot_ZJetBalance() {

  if(numFiles != nGen) 
    cout << "Not All Files Are Connected to Chain: Something Wrong !" << endl;

  cout<< "nentries = " << mychain->GetEntries() << endl;

  //  return;


  // ****** get the weight value for each pthat bin ******* //
  for(int i=0; i<nGen; i++) {
    TFile f( basedir + inFileNames[i], "read");
    TTree* atree = (TTree*) f.Get("ZJet");
    nEvents[i]  = (int) atree->GetEntries();
    if (applyWeights) weight[i] = lumi * crosssection[i]  /  nEvents[i];
    else weight[i] = 1.0;
    cout << "weight_" << i << " == " << weight[i] << endl;
  }

  for(int i=0; i<nZPtBins; i++) {
    pt[i] = 0.5*( theZPt[i] + theZPt[i+1] );
    ptBin[i] = Form("%d_%d", (int) theZPt[i], (int) theZPt[i+1] );
    errpt[i]     = 0.0;
  }

  // *************** get all the data branches *************** //
  SetupTreeBranches( mychain );


//   for(int i=0; i<5; i++) {plot_ZJetBalance(i);}
  for(int i=1; i<5; i++) {plot_ZJetBalance(i);}
}









////////////////// pT balance for a given jet algo ////////

void plot_ZJetBalance(int JetAlgo) {


  Float_t genMean[nZPtBins];
  Float_t genSigma[nZPtBins]; 
  Float_t recoMean[nZPtBins]; 
  Float_t recoSigma[nZPtBins];

  std::string algoSt;
  if(JetAlgo==0) algoSt = "_ic5";
  if(JetAlgo==1) algoSt = "_sc5";
  if(JetAlgo==2) algoSt = "_sc7";
  if(JetAlgo==3) algoSt = "_kt4";
  if(JetAlgo==4) algoSt = "_kt6";

  std::string textFilename = "ZjetPtBalance" + algoSt + ".txt";
  std::string JetResponseFilename 
    = "Histograms_ZjetResponse" + algoSt + ".root";


  
  plot_ZJetBalance(  JetAlgo, JetResponseFilename, genMean, 
		     genSigma, recoMean, recoSigma );


  //////////////////////////////////////////////////
  if(storeResponseInTextFile) {
    
    FILE *file = fopen(textFilename.c_str(),"w+");
    fprintf( file ,"%s    %s   %s   %s   %s \n", "pT bin", 
	     "genMean", "genSigma", "recoMean", "recoSigma");
    
    for(int i=0; i<nZPtBins; i++) {
      fprintf( file ,"%5.1f   %5.4f   %5.4f   %5.4f   %5.4f \n", pt[i],
	       genMean[i], genSigma[i], recoMean[i], recoSigma[i] );
    }
    fclose(file); 
  }
  //////////////////////////////////////////////////

  if(makeplot_Response) 
    DrawResponsePlot( genMean, genSigma, recoMean, recoSigma);
}










////////////////// pT balance for a given Z pT bin ////////


void plot_ZJetBalance( int index, std::string JetResponseFilename, 
		       Float_t genMean[], Float_t genSigma[], 
		       Float_t recoMean[], Float_t recoSigma[]) {


  ////////////// Defining the L3JetCorrector ///////////////////////
  double p[6];
  p[0] = 10.0;          
  p[1] = 800.0;         

  // // with MEAN response
//   p[2] = 0.953124;
//   p[3] = 4.87151;
//   p[4] = 2.83723;  
//   p[5] = 2.91468;


  // // with MPV response
//   p[2] = 0.978755;
//   p[3] = 2.1759;
//   p[4] = 2.25031;  
//   p[5] = 0.0;


  // // with TruncMEAN 2 sigma response
//   p[2] = 1.00299;
//   p[3] = 3.83695;
//   p[4] = 2.87351;  
//   p[5] = 1.6071;


  // // with TruncMEAN 1.5 sigma response
  p[2] = 1.00299;
  p[3] = 3.83695;
  p[4] = 2.87351;  
  p[5] = 1.6071;


  // // with TruncMEAN 1.5 sigma iterative 3 times
//   p[2] = 1.16793;
//   p[3] = 8.8575;
//   p[4] = 5.01816;  
//   p[5] = 5.03132;


  // // with TruncMEAN 1.5 sigma iterative n times
//   p[2] = 1.23392;
//   p[3] = 9.33578;
//   p[4] = 5.60468;  
//   p[5] = 5.52937;


  // // with TruncMEAN 1 sigma response
//   p[2] = 0.690152;
//   p[3] = 2.97111;
//   p[4] = 1.70586;  
//   p[5] = 0.72181;


  // // with MPV: when p[5] is allowed to go negative       
//   p[2] = 0.893141;
//   p[3] = 0.00355417;
//   p[4] = 0.0123656;  
//   p[5] = -1.00176;


  // // from dijet MC truth      
//   p[2] = 0.996998;
//   p[3] = 4.39412;
//   p[4] = 2.96134;  
//   p[5] = 1.69966;


   if(storeResponsHistogramsInRootFile == true) {
     TFile respHistFile(JetResponseFilename.c_str(),"RECREATE");
  }

  TH1F* responseHistGen[nZPtBins];
  TH1F* responseHistReco[nZPtBins];
  TH1F* genJetpt[nZPtBins];
  TH1F* recoZpt[nZPtBins];
  TH1F* caloJetpt[nZPtBins];
  TH1F* Zpt[nZPtBins];

  for(int i=0; i<nZPtBins; i++) {

    responseHistGen[i] = new TH1F(TString("responseHistGen_")+ptBin[i],
				  "",40, 0.0, 2.0);
    TAxis* responseHistGenx = responseHistGen[i]->GetXaxis();
    TAxis* responseHistGeny = responseHistGen[i]->GetYaxis();
    responseHistGenx->SetTitle("p_{T}^{jet} / p_{T}^{Z}   ");
    responseHistGeny->SetTitle("Events / 0.05");
    responseHistGeny->SetTitleOffset(1.4);   
    responseHistGenx->SetNdivisions(505);
    responseHistGeny->SetNdivisions(505);

    responseHistReco[i] = new TH1F(TString("responseHistReco_")+ptBin[i],
				   "",40,0.0,2.0);
    responseHistReco[i]->SetLineColor(2);
    responseHistReco[i]->SetMarkerColor(2);    
    genJetpt[i] = new TH1F(TString("genJetpt_")+ptBin[i],
			   "", 140, 0, 800);
    caloJetpt[i] = new TH1F(TString("caloJetpt_")+ptBin[i],
			    "", 140, 0, 800);
    recoZpt[i] = new TH1F(TString("recoZpt_")+ptBin[i],
			  "", 140, 0, 800);

    responseHistReco[i]->Sumw2();
    responseHistGen[i]->Sumw2();
    genJetpt[i]->Sumw2();
    caloJetpt[i]->Sumw2();
    recoZpt[i]->Sumw2();
  }


  for (Long64_t entry =0; entry < mychain->GetEntries(); entry++) {
    
    mychain->GetEntry(entry);
    if(entry%100000==0) std::cout<<"**** Event # "<< entry <<std::endl;


    // Fill generator level quantities
    int leadGenIndex=-1, secondGenIndex=-1;
    FindLeadIndex(JetGenPt[index], JetGenEta[index], JetGenPhi[index],
		  eMinusEtaGen, eMinusPhiGen, 
		  ePlusEtaGen, ePlusPhiGen, leadGenIndex, secondGenIndex);

    if(leadGenIndex !=-1 && secondGenIndex !=-1 && 
       (ePlusPtGen>20.0) && (eMinusPtGen>20.0) && 
       ((fabs(ePlusEtaGen)<1.4442) || 
	(fabs(ePlusEtaGen)>1.560 && fabs(ePlusEtaGen)<2.5)) && 
       ((fabs(eMinusEtaGen)<1.4442) || 
	(fabs(eMinusEtaGen)>1.560 && fabs(eMinusEtaGen)<2.5)) && 
       (fabs(JetGenEta[index][leadGenIndex])<1.3) ) {
	
      float leadGenJetPt   = JetGenPt[index][leadGenIndex];
      float secondGenJetPt = JetGenPt[index][secondGenIndex];

      double wt = GetWeight( mychain->GetFile()->GetName() );

      double ptRatioGen  = (double) leadGenJetPt/ (double) Z_PtGen;
      if((double) secondGenJetPt/(double) Z_PtGen < 0.1) {  

	for(int bin=0; bin<nZPtBins; bin++) { //begin Z pT bin loop
	  if(  (Z_PtGen > theZPt[bin]) && (Z_PtGen < theZPt[bin+1]) ) {
	    responseHistGen[bin]->Fill( ptRatioGen, wt );
	    genJetpt[bin]->Fill(leadGenJetPt, wt );  
	  } // end if loop
	} // end Z pT loop
      }
    }


    // Fill reco level quantities
    int leadRecoIndex=-1, secondRecoIndex=-1;
    FindLeadIndex(JetRecoPt[index], JetRecoEta[index], JetRecoPhi[index],
		  eMinusEta, eMinusPhi, 
		  ePlusEta, ePlusPhi, leadRecoIndex, secondRecoIndex);
    
    if(leadRecoIndex !=-1 && secondRecoIndex !=-1) {

      float leadRecoJetPt;
      float secondRecoJetPt;
      float leadRecoJetEta;
      double dPhiReco;

      double L3scale = 1.0;
      double apt = JetRecoPt[index][leadRecoIndex];
      if( apt< p[0] ) apt = p[0];
      if( apt> p[1] ) apt = p[1];
      double log10pt = log10(apt);
      double result = p[2]+p[3]/(pow(log10pt,p[4])+p[5]);

      if(usingCorrectedCaloJetPt) L3scale = result;

      leadRecoJetPt   = L3scale * JetRecoPt[index][leadRecoIndex];
      secondRecoJetPt = JetRecoPt[index][secondRecoIndex];
      leadRecoJetEta  = JetRecoEta[index][leadRecoIndex];
      dPhiReco = dPhi(JetRecoPhi[index][leadRecoIndex], Z_Phi);


      bool pass = BoolCutResult( Z_Pt, mZee, eMinusPt, eMinusEta, ePlusPt, 
				 ePlusEta, eMinus_trackiso, ePlus_trackiso, 
				 eMinus_ecaliso, ePlus_ecaliso, 
				 eMinus_hcaliso, ePlus_hcaliso, 
				 iseMinusLoose, isePlusLoose,
				 leadRecoJetPt, leadRecoJetEta, 
				 secondRecoJetPt, dPhiReco);

      if( ! pass )  continue;

      double wt = GetWeight( mychain->GetFile()->GetName() );
      double ptRatioReco = (double) leadRecoJetPt/ (double) Z_Pt;
 

      for(int bin=0; bin<nZPtBins; bin++) { //begin Z pT bin loop
	if(  (Z_Pt > theZPt[bin]) && (Z_Pt < theZPt[bin+1]) ) {
	  responseHistReco[bin]->Fill( ptRatioReco, wt  );
	  caloJetpt[bin]->Fill(leadRecoJetPt, wt );
	  recoZpt[bin]->Fill( Z_Pt, wt  );   
	} // end if loop
      } // end Z pT loop
    }
  } // end TTree loop


  // Fill the mean and error vectors

    for(int i=0; i<nZPtBins; i++) {
      if(ScaleToLuminosity) {
	for(int j=0; j<40; j++) {
	  double err = sqrt( responseHistReco[i]->GetBinContent(j));
	  responseHistReco[i]->SetBinError(j, err);
	  err = sqrt( responseHistGen[i]->GetBinContent(j));
	  responseHistGen[i]->SetBinError(j, err);
	}
      }

      genMean[i]   = (Float_t) responseHistGen[i]->GetMean(1);
      genSigma[i]  = (Float_t) responseHistGen[i]->GetMean(11);
      recoMean[i]  = (Float_t) responseHistReco[i]->GetMean(1); 
      recoSigma[i] = (Float_t) responseHistReco[i]->GetMean(11);
    }



  // plot the pT balance response histograms
  if(makeplot_ZptBalance==true) {

    for(int i=0; i<nZPtBins; i++) {
      TString plotname = Form("ptBalance-allCuts_%d_%d", (int) theZPt[i], 
			      (int) theZPt[i+1] );  
      PlotOnCanvas( *responseHistGen[i], *responseHistReco[i], plotname);
    }
  }   


  // Now write all the histograms in a ROOT file
  if(storeResponsHistogramsInRootFile == true) {
    respHistFile.cd();
    for(int i=0; i<nZPtBins; i++) {
      responseHistGen[i]->Write();
      responseHistReco[i]->Write();
      genJetpt[i]->Write();
      caloJetpt[i]->Write();
      recoZpt[i]->Write();
    }
    respHistFile.Close();
  }


  // clean up the memory
  delete [] responseHistGen;
  delete [] responseHistReco;
  delete [] genJetpt;
  delete [] caloJetpt;
  delete [] recoZpt;
}






////////// Apply event selection cuts ///////////////////

bool BoolCutResult( float zPt, float mZ, float e1Pt, 
		    float e1Eta, float e2Pt, float e2Eta, 
		    float e1trackiso, float e2trackiso, 
		    float e1ecaliso, float e2ecaliso, 
		    float e1hcaliso, float e2hcaliso, 
		    bool ise1Loose, bool ise2Loose,
		    float leadPt, float leadEta, float secondPt, 
		    float phidiff) {

  bool result = true;

  // Z mass cut
  if( fabs(mZ-91.2) > 10 ) result = false;


  // electron pT cut
  if( e1Pt < 20.0 ) result = false;
  if( e2Pt < 20.0 ) result = false;


  // electron acceptance
  if( !((fabs(e1Eta)<1.4442) || 
	(fabs(e1Eta)>1.560 && fabs(e1Eta)<2.5)) ) result = false;
  if( !((fabs(e2Eta)<1.4442) || 
	(fabs(e2Eta)>1.560 && fabs(e2Eta)<2.5)) ) result = false;


  // electron isolation: e1 barrel
  if( fabs(e1Eta)<1.4442) { 
    if(e1trackiso > 7.2 )  result = false; 
    if(e1ecaliso > 5.7 )  result = false; 
    //   if(e1hcaliso > 8.1 )  result = false; 
  }

  // e1 endcap
  if( fabs(e1Eta)>1.560 && fabs(e1Eta)<2.5 ) { 
    if(e1trackiso > 5.1 )  result = false; 
    if(e1ecaliso > 5.0 )  result = false; 
    // if(e1hcaliso > 3.4 )  result = false; 
    }

  // e2 barrel
  if( fabs(e2Eta)<1.4442) { 
    if(e2trackiso > 7.2 )  result = false; 
    if(e2ecaliso > 5.7 )  result = false;
    // if(e2hcaliso > 8.1 )  result = false; 
  }

  // e2 endcap
  if( fabs(e2Eta)>1.560 && fabs(e2Eta)<2.5 ) { 
      if(e2trackiso > 5.1 )  result = false; 
      if(e2ecaliso > 5.0 )  result = false; 
      //   if(e2hcaliso > 3.4 )  result = false; 
  }

  // electron Id

//   if(!(ise1Loose)) result = false; 
//   if(!(ise2Loose)) result = false; 

  // Minimum jet pT cut on all jets
  if( leadPt < MINPTCUT ) result = false;
  if( secondPt < MINPTCUT ) result = false;

  // Cut on the second jet pT
  if( secondPt / zPt > 0.1 ) result = false;

  // Leading jet has to be in the barrel
  if( fabs(leadEta) > 1.3 ) result = false;

  // Z and lead jet are back-to-back in phi
  if( fabs(phidiff) < 2.94 ) result = false;

  return result;
}







// Draw pT balance plots
void PlotOnCanvas(TH1& genHist, TH1& recoHist, TString plotname) {

  gROOT->ProcessLine(".L mystyle.C");
  setTDRStyle();
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPadLeftMargin(0.14);
  tdrStyle->SetPadRightMargin(0.10);
  tdrStyle->SetLegendBorderSize(0);
  tdrStyle->SetTitleYOffset(1.5);

  TCanvas canvas("canvas", "", 500, 500);
  gStyle->SetOptStat(0);
  genHist.SetMinimum(0);
  genHist.Draw();
  genHist.Draw("hist same");
  recoHist.Draw("same");
  recoHist.Draw("HIST same");
  leg_hist = new TLegend(0.6,0.7,0.89,0.89);
  leg_hist->AddEntry(&genHist,"Generator level","l");
  leg_hist->AddEntry(&recoHist,"Calorimeter level","l");
  leg_hist->SetFillColor(0);
  leg_hist->Draw();

  canvas.SaveAs(plotname+TString(".eps"));
  canvas.SaveAs(plotname+TString(".gif"));
  canvas.SaveAs(plotname+TString(".root"));
  canvas.Close();

//   delete leg_hist;
}






// ******************************************************* //
// ********* compare Zee response with Zmumu *************** //


void DrawResponsePlot( Float_t genMean[], Float_t genSigma[], 
		       Float_t recoMean[], Float_t recoSigma[]) {


  gROOT->ProcessLine(".L mystyle.C");
  setTDRStyle();
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPadLeftMargin(0.2);
  tdrStyle->SetPadRightMargin(0.10);
  tdrStyle->SetLegendBorderSize(0);
  tdrStyle->SetTitleYOffset(1.3);

  // plot full spectrum
  TGraphErrors *ptbalanceGen  = new TGraphErrors(nZPtBins, pt, genMean, 
						 errpt, genSigma);
  TGraphErrors *ptbalanceReco = new TGraphErrors(nZPtBins, pt, recoMean, 
						 errpt, recoSigma);
  // plot Zmumu values
  Float_t ptmm[9] = { 40.0, 60.0, 100.0, 140.0, 200.0, 250.0, 
		      330.0, 400.0, 520.0 };
  Float_t balancemm[9] = { 0.496, 0.568, 0.66, 0.71, 0.75, 0.765, 
			   0.775, 0.79, 0.81 }; 
  TGraph *ptbalancemm = new TGraph( 9, ptmm, balancemm);


  ptbalanceGen->GetXaxis()->SetTitle("p_{T}^{Z} (GeV/c)");
  ptbalanceGen->GetYaxis()->SetTitle("p_{T}^{jet} / p_{T}^{Z}");
  ptbalanceGen->SetMarkerStyle(22);
  ptbalanceGen->SetMarkerSize(1.2);
  ptbalanceGen->SetTitle("");
  ptbalanceGen->SetMinimum(0.3);
  ptbalanceReco->SetMarkerColor(2);
  ptbalanceReco->SetLineColor(2);
  ptbalanceReco->SetMarkerStyle(22);
  ptbalanceReco->SetMarkerSize(1.2);
  ptbalanceReco->SetMinimum(0.3);
  ptbalancemm->SetMarkerStyle(24);
  ptbalancemm->SetMarkerSize(1.2);
  ptbalancemm->SetMinimum(0.3);
  ptbalancemm->SetMarkerColor(4);
  ptbalancemm->SetLineColor(4);

  TCanvas c1("c1","",500,500);
  ptbalanceGen->Draw("APL");
  ptbalanceReco->Draw("PL");
  ptbalancemm->Draw("PL");
  leg_hist = new TLegend(0.6,0.7,0.89,0.89);
  leg_hist->AddEntry( ptbalanceGen, "Generator level", "l");
  leg_hist->AddEntry( ptbalanceReco,"Calorimeter level","l");
  leg_hist->AddEntry( ptbalancemm,"Z#rightarrow#mu#mu","l");
  leg_hist->SetFillColor(0);
  leg_hist->Draw();
  c1.SaveAs("PtBalanceVsPt.eps");
  c1.SaveAs("PtBalanceVsPt.gif");
  c1.SaveAs("PtBalanceVsPt.root");
  c1.Close();
  delete leg_hist;
  delete ptbalanceGen;
  delete ptbalanceReco;
  delete ptbalancemm;
}



// ******************************************************* //
// ******************************************************* //
// ******************************************************* //
// ******************************************************* //

void SetupTreeBranches( TTree* tree) {

  tree->SetBranchAddress("JetGenPt",   JetGenPt);
  tree->SetBranchAddress("JetGenEta",  JetGenEta);
  tree->SetBranchAddress("JetGenPhi",  JetGenPhi);
  tree->SetBranchAddress("JetCorPt",   JetCorPt);
  tree->SetBranchAddress("JetCorEta",  JetCorEta);
  tree->SetBranchAddress("JetCorPhi",  JetCorPhi);
  tree->SetBranchAddress("JetRecoPt",  JetRecoPt);
  tree->SetBranchAddress("JetRecoEta", JetRecoEta);
  tree->SetBranchAddress("JetRecoPhi", JetRecoPhi);

  tree->SetBranchAddress("Z_Pt",       &Z_Pt);
  tree->SetBranchAddress("Z_Phi",      &Z_Phi);
  tree->SetBranchAddress("Z_Eta",      &Z_Eta);
  tree->SetBranchAddress("Z_PtGen",    &Z_PtGen);
  tree->SetBranchAddress("Z_PhiGen",   &Z_PhiGen);
  tree->SetBranchAddress("ePlusPt",    &ePlusPt);
  tree->SetBranchAddress("eMinusPt",   &eMinusPt);
  tree->SetBranchAddress("ePlusPtGen", &ePlusPtGen);
  tree->SetBranchAddress("eMinusPtGen",&eMinusPtGen);
  tree->SetBranchAddress("ePlusEta",   &ePlusEta);
  tree->SetBranchAddress("eMinusEta",  &eMinusEta);
  tree->SetBranchAddress("ePlusEtaGen",   &ePlusEtaGen);
  tree->SetBranchAddress("eMinusEtaGen",  &eMinusEtaGen);
  tree->SetBranchAddress("ePlusPhi",   &ePlusPhi);
  tree->SetBranchAddress("eMinusPhi",  &eMinusPhi);
  tree->SetBranchAddress("ePlusPhiGen",   &ePlusPhiGen);
  tree->SetBranchAddress("eMinusPhiGen",  &eMinusPhiGen);
  tree->SetBranchAddress("mZeeGen",     &mZeeGen);
  tree->SetBranchAddress("mZee",        &mZee);
  tree->SetBranchAddress("ePlusE",            &ePlusE);
  tree->SetBranchAddress("eMinusE",           &eMinusE);
  tree->SetBranchAddress("eMinus_ecaliso",     &eMinus_ecaliso);
  tree->SetBranchAddress("ePlus_ecaliso",      &ePlus_ecaliso);
  tree->SetBranchAddress("eMinus_hcaliso",     &eMinus_hcaliso);
  tree->SetBranchAddress("ePlus_hcaliso",      &ePlus_hcaliso);
  tree->SetBranchAddress("ePlus_trackiso",  &ePlus_trackiso);	   
  tree->SetBranchAddress("eMinus_trackiso", &eMinus_trackiso);
  tree->SetBranchAddress("iseMinusLoose", &iseMinusLoose);
  tree->SetBranchAddress("isePlusLoose", &isePlusLoose);
  tree->SetBranchAddress("NumRecoJets",  &NumRecoJets);
  tree->SetBranchAddress("NumGenJets",   &NumGenJets);  

  tree->SetBranchStatus("*",    0);
  tree->SetBranchStatus("JetGenPt",    1);
  tree->SetBranchStatus("JetGenEta",   1);
  tree->SetBranchStatus("JetGenPhi",   1);
  tree->SetBranchStatus("JetRecoPt",   1);
  tree->SetBranchStatus("JetRecoEta",  1);
  tree->SetBranchStatus("JetRecoPhi",  1);
  tree->SetBranchStatus("JetCorPt",    1);
  tree->SetBranchStatus("JetCorEta",   1);
  tree->SetBranchStatus("JetCorPhi",   1);
  tree->SetBranchStatus("Z_Pt",        1);
  tree->SetBranchStatus("Z_Phi",       1);
  tree->SetBranchStatus("Z_Eta",       1);
  tree->SetBranchStatus("Z_PtGen",     1);
  tree->SetBranchStatus("Z_PhiGen",    1);
  tree->SetBranchStatus("ePlusPt",     1);
  tree->SetBranchStatus("eMinusPt",    1);
  tree->SetBranchStatus("ePlusPtGen",  1);
  tree->SetBranchStatus("eMinusPtGen", 1);
  tree->SetBranchStatus("ePlusEta",    1);
  tree->SetBranchStatus("eMinusEta",   1);
  tree->SetBranchStatus("ePlusEtaGen",    1);
  tree->SetBranchStatus("eMinusEtaGen",   1);
  tree->SetBranchStatus("ePlusPhi",    1);
  tree->SetBranchStatus("eMinusPhi",   1);
  tree->SetBranchStatus("ePlusPhiGen",    1);
  tree->SetBranchStatus("eMinusPhiGen",   1);
  tree->SetBranchStatus("mZeeGen",     1);
  tree->SetBranchStatus("mZee",        1);
  tree->SetBranchStatus("ePlusE",            1);
  tree->SetBranchStatus("eMinusE",           1);
  tree->SetBranchStatus("ePlus_ecaliso",      1);
  tree->SetBranchStatus("eMinus_ecaliso",     1);
  tree->SetBranchStatus("ePlus_hcaliso",      1);
  tree->SetBranchStatus("eMinus_hcaliso",     1);
  tree->SetBranchStatus("ePlus_trackiso",  1);	   
  tree->SetBranchStatus("eMinus_trackiso", 1);
  tree->SetBranchStatus("iseMinusLoose", 1);
  tree->SetBranchStatus("isePlusLoose",  1);
  tree->SetBranchStatus("NumRecoJets",  1);
  tree->SetBranchStatus("NumGenJets",  1);  
}





// ******************************************************* //
// ******************************************************* //
// ******************************************************* //
// ******************************************************* //


double GetWeight(TString filename){

  double wt =0.0;
  for(int i=0; i<nGen; i++) {
    TString st2 = basedir + inFileNames[i];
    if(  filename.CompareTo(st2) == 0) wt = weight[i]; 
  }

  return wt;
}
// ******************************************************* //




// find the leading and second jet indices

void FindLeadIndex(float pT[], float eta[], float phi[],
		   float e1eta, float e1phi, 
		   float e2eta, float e2phi, int &lead, int &sec) {
  float max = 0.0;
  lead = -1;
  for (int i=0; i<10; i++) {
    
    double r1 = radius( eta[i],phi[i], e1eta, e1phi );
    double r2 = radius( eta[i],phi[i], e2eta, e2phi );
    if(!(r1>0.2 && r2>0.2) ) continue;
    if(pT[i]>max) { max = pT[i]; lead = i; }
  }

  max = 0.0;
  sec = -1;
  for (int i=0; i<10; i++) {
    if(i==lead) continue;
    double r1 = radius( eta[i],phi[i], e1eta, e1phi );
    double r2 = radius( eta[i],phi[i], e2eta, e2phi );
    if(!(r1>0.2 && r2>0.2) ) continue;
    if(pT[i]>max) { max = pT[i]; sec = i; }
  }
}




void FindLeadIndices( int algo, int &leadReco, int &secReco, 
		      int &leadGen, int &secGen ) {

  FindLeadIndexReco( algo, leadReco, secReco );
  FindLeadIndexGen( algo, leadGen, secGen );
} 



void FindLeadIndexReco( int algo, int &lead, int &sec) {

  float eta, phi;

  // calculate dR of all the jets w.r.t. the two electrons
  for( int j=0; j < NumRecoJets; j++ ) {

    if(usingCorrectedCaloJetPt) {
      p[j]  = JetCorPt[algo][j];
      eta = JetCorEta[algo][j];
      phi = JetCorPhi[algo][j];
    }
    else {
      p[j]  = JetRecoPt[algo][j];
      eta = JetRecoEta[algo][j];
      phi = JetRecoPhi[algo][j];
    }

    // protection:
    if( p[j] < 0.001 ) { p[j] = 0.001; break; }

    R1[j] = radius( eta, phi, eMinusEta,  eMinusPhi );
    R2[j] = radius( eta, phi, ePlusEta,  ePlusPhi );
   }
  
  FindIndex( p, R1, R2, lead, sec, NumRecoJets );
}





void FindLeadIndexGen( int algo, int &lead, int &sec ) {

  for( int j=0; j < NumGenJets; j++ ) {

    p[j]  = JetGenPt[algo][j];

    // protection:
    if( p[j] < 0.001 ) { p[j] = 0.001; break; }

    R1[j] = radius( JetGenEta[algo][j], JetGenPhi[algo][j], 
			eMinusEtaGen,  eMinusPhiGen );
    R2[j] = radius( JetGenEta[algo][j], JetGenPhi[algo][j], 
			ePlusEtaGen,  ePlusPhiGen );
  }
  
  FindIndex( p, R1, R2, lead, sec, NumGenJets );
}





void FindIndex( float pT[], float dr1[], float dr2[], int &lead, 
		int &sec, int maxIndx = NUM_JET_MAX ) {
  float max = 0.0;
  lead = -1; 
  sec = -1;

  for (int i=0; i< maxIndx; i++) {
    if( !(dr1[i]>0.1 && dr2[i]>0.1) ) continue;
    if( pT[i]>max ) {  max = pT[i];  lead = i; }
  }

  if( lead == -1 ) return;

  max = 0.0;
  for (int i=0; i< maxIndx; i++) {
    if( !(dr1[i]>0.1 && dr2[i]>0.1) || i==lead) continue;
    if( pT[i]>max ) {  max = pT[i];  sec = i; }
  }

}




