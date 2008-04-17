
{
const int NPtBins=11;
const int NETA = 1;
 const double Pt[NPtBins+1] = {20,25,30,40,50,60,70,100,140,200,400,800};
const double eta_boundaries[NETA+1] = {-1.3,1.3};

  TFile *f;
  TH1F *hResponse,*hMeanRefPt,*hMeanCaloPt;
  double yRefPt[NPtBins],eyRefPt[NPtBins],xRefPt[NPtBins],exRefPt[NPtBins];
  double yCaloPt[NPtBins],eyCaloPt[NPtBins],xCaloPt[NPtBins],exCaloPt[NPtBins];
  double x,y,ex,ey,x1,ex1;
  int i,N;
  f = new TFile("FitterResults_test_Icone5.root","r");
  if (f->IsZombie()) break; 
  hResponse = (TH1F*)f->Get("Response");
  hMeanRefPt = (TH1F*)f->Get("MeanRefPt");
  hMeanCaloPt = (TH1F*)f->Get("MeanCaloPt");
  N = 0;
  for(i=0;i<NPtBins;i++)
    {
      y = hResponse->GetBinContent(i+1);
      ey = hResponse->GetBinError(i+1);
      x = hMeanRefPt->GetBinContent(i+1);
      ex = hMeanRefPt->GetBinError(i+1);
      x1 = hMeanCaloPt->GetBinContent(i+1);
      ex1 = hMeanCaloPt->GetBinError(i+1); 
      if (y>0 && x>0 && x1>0 && ey>0.000001 && ey<0.2)
	{
          yRefPt[N] = y;
          eyRefPt[N] = ey;
          xRefPt[N] = x;
          exRefPt[N] = ex;
          xCaloPt[N] = x1;
          exCaloPt[N] = ex1;
	  N++;
	}  
    }
  TGraphErrors *ptbalanceReco = new TGraphErrors(N,xRefPt,yRefPt,exRefPt,eyRefPt);
  ptbalanceReco->SetMarkerStyle(20);
  ptbalanceReco->SetLineColor(2);
  ptbalanceReco->SetMarkerColor(2);

  gROOT->ProcessLine(".L mystyle.C");
  setTDRStyle();
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPadLeftMargin(0.2);
  tdrStyle->SetPadRightMargin(0.10);
  tdrStyle->SetLegendBorderSize(0);
  tdrStyle->SetTitleYOffset(1.3);

  // plot full spectrum
  //  TGraphErrors *ptbalanceReco = new TGraphErrors(NGenPtBins, pt, recoMean, 
  //						 errpt, recoSigma);
  // plot Zmumu values
  Float_t ptmm[9] = { 40.0, 60.0, 100.0, 140.0, 200.0, 250.0, 
		      330.0, 400.0, 520.0 };
  Float_t balancemm[9] = { 0.496, 0.568, 0.66, 0.71, 0.75, 0.765, 
			   0.775, 0.79, 0.81 }; 
  TGraph *ptbalancemm = new TGraph( 9, ptmm, balancemm);
  ptbalancemm->SetMarkerColor(4);
  ptbalancemm->SetLineColor(4);


//   Float_t ptee07[9] = { 45.0, 70.0, 100.0, 145.0, 200.0, 265.0, 
// 			340.0, 425.0, 535.0}; 
//   Float_t ptee07Err[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

//   Float_t balance07ee[9] = { 0.6749, 0.7079, 0.7120, 0.7453,
// 			     0.7789, 0.8014, 0.8162, 0.8184, 0.8286 }; 
//   Float_t balance07eeErr[9] = { 0.0266, 0.0208, 0.0138, 0.0087, 
// 				0.0063, 0.0053, 0.0044, 0.0033, 0.0028 }; 


//   Float_t ptee07[8] = { 34.2216, 59.6234, 82.7739, 116.84, 164.222, 267.102, 
// 			425.0, 535.0}; 
//   Float_t ptee07Err[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

//   Float_t balance07ee[8] = { 0.497503, 0.596338, 0.627716, 0.700546,
// 			     0.740916, 0.779794, 0.8084, 0.8186 }; 
//   Float_t balance07eeErr[8] = { 0.0175227, 0.0220648, 0.0230714,
// 				0.00509754, 0.00352788, 0.00190956,
// 				0.0033, 0.0028 }; 
//   TGraphErrors *ptbalance07ee = new TGraphErrors(8, ptee07, balance07ee, 
// 						 ptee07Err, balance07eeErr);
//   ptbalance07ee->SetMarkerColor(2);
//   ptbalance07ee->SetLineColor(2);


  // TFile fZmmJ08("corrections_responses_Summer08_Zmumujet_100pb.root","read");
 
  TFile fZmmJ08( "ZPJ_responses_vs_PtZ_100pb.root" );
//  TGraphErrors *ptbalance08mm = 
//     (TGraphErrors *) fZmmJ08.Get("hm_L2CorJetIC5Calo_ratio_graph"); 
//  TGraphErrors *ptbalance08mm = 
//     (TGraphErrors *) fZmmJ08.Get("L2CorJetIC5Calo"); 




 TGraphErrors *ptbalance08mm = 
    (TGraphErrors *) fZmmJ08.Get("iterativeCone5CaloJets"); 
  ptbalance08mm->SetMarkerColor(7);
  ptbalance08mm->SetLineColor(7);
 

//   TFile fGammaJ07("csa07.root","read");
//   TGraphErrors *ptbalance07Gamma = 
//     (TGraphErrors *) fGammaJ07.Get("response_ite"); 
//   ptbalance07Gamma->SetMarkerColor(6);
//   ptbalance07Gamma->SetLineColor(6);
//   ptbalance07Gamma->SetMarkerSize(0.8);


//   TFile fGammaJ08("csa08.root","read");
//   TGraphErrors *ptbalance08Gamma = 
//     (TGraphErrors *) fGammaJ08.Get("response_ite"); 
//   ptbalance08Gamma->SetMarkerColor(8);
//   ptbalance08Gamma->SetLineColor(8);
//   ptbalance08Gamma->SetMarkerSize(0.9);


  TFile fGammaJSumm08("summer08.root","read"); 
//   TGraphErrors *ptbalanceSumm08Gamma = 
//     (TGraphErrors*)fGammaJSumm08->Get("resp_sig");

//   TGraphErrors *ptbalanceSumm08Gamma = 
//     (TGraphErrors*)fGammaJSumm08->Get("measrespvsptphot_sig");


  TGraphErrors *gammaJetResp = 
    (TGraphErrors*)fGammaJSumm08.Get("measrespvsptphot_sig");
  TGraphErrors *gammaJetCorr = 
    (TGraphErrors*)fGammaJSumm08.Get("corr_sig_100pb");

  Double_t* vectX = gammaJetResp->GetX(); 
  Double_t* vectEX = gammaJetResp->GetEX();
  Double_t* vectY = gammaJetResp->GetY(); 
  Double_t* vectEY = gammaJetCorr->GetEY();
  int mumPts = gammaJetCorr->GetN();

  for (int j=0; j<mumPts; ++j) {
    vectEY[j] = vectEY[j] / pow(vectY[j], 2); 
  }
  TGraphErrors *ptbalanceSumm08Gamma 
    = new TGraphErrors(mumPts,vectX,vectY,vectEX,vectEY);



  ptbalanceSumm08Gamma->SetMarkerColor(1);
  ptbalanceSumm08Gamma->SetLineColor(1);
  ptbalanceSumm08Gamma->SetMarkerStyle(20);
  //  TF1 *cc1 = (TF1*)ptbalanceSumm08Gamma->GetFunction("fsig_b");
  TF1 *cc1 = (TF1*)ptbalanceSumm08Gamma->GetFunction("fsig");
  delete cc1;
  TPaveStats *p1 
    = (TPaveStats*) ptbalanceSumm08Gamma->GetListOfFunctions()->FindObject("stats");
  ptbalanceSumm08Gamma->GetListOfFunctions()->Remove(p1);

  TAxis* xaxis = ptbalanceSumm08Gamma->GetXaxis();
  TAxis* yaxis = ptbalanceSumm08Gamma->GetYaxis();
  xaxis->SetTitle("p_{T}^{Z/#gamma} (GeV/c)");
  yaxis->SetTitle("p_{T}^{jet} / p_{T}^{Z/#gamma}   ");
  xaxis->SetMoreLogLabels();
  xaxis->SetNoExponent();
  ptbalanceSumm08Gamma->SetTitle("");
  ptbalanceSumm08Gamma->SetMinimum(0.3);
  ptbalanceSumm08Gamma->SetMaximum(1);

  ptbalanceReco->RemovePoint(9);
  ptbalance08mm->RemovePoint(0);


  TCanvas c1("c1","",500,500);
  ptbalanceSumm08Gamma->Draw("APL");
  ptbalanceReco->Draw("PL");
  // ptbalance07ee->Draw("PL");
  // ptbalance07Gamma->Draw("PL");
  // ptbalance08Gamma->Draw("PL");

  // ptbalancemm->Draw("PL");
  ptbalance08mm->Draw("PL");
  leg_hist = new TLegend(0.38,0.25,0.88,0.5);
  leg_hist->AddEntry( ptbalanceReco,"Summer08 Z#rightarrow ee + jet","P");
  //  leg_hist->AddEntry( ptbalance07ee,"CSA07 Z#rightarrow ee + jet","P");
  //  leg_hist->AddEntry( ptbalancemm,"CSA07 Z#rightarrow#mu#mu + jet","P");
  leg_hist->AddEntry( ptbalance08mm,"Summer08 Z#rightarrow#mu#mu + jet","P");
  // leg_hist->AddEntry( ptbalance07Gamma,"CSA07 #gamma + jet","P");
  // leg_hist->AddEntry( ptbalance08Gamma,"CSA08 #gamma + jet","P");
  leg_hist->AddEntry( ptbalanceSumm08Gamma,"Summer08 #gamma + jet","P");
  leg_hist->SetMargin(0.1);
  leg_hist->SetFillColor(0);
  leg_hist->Draw();
  gPad->SetLogx();
  c1.SaveAs("PtBalanceVsPt.eps");
  c1.SaveAs("PtBalanceVsPt.gif");
  c1.SaveAs("PtBalanceVsPt.root");
  //   c1.Close();
  //   delete leg_hist;
  //   delete ptbalanceGen;
  //   delete ptbalanceReco;
  //   delete ptbalancemm;

}
