void plotClosure() {

  gROOT->ProcessLine(".L mystyle.C");
  setTDRStyle();
  tdrStyle->SetPadLeftMargin(0.18);
  tdrStyle->SetPadRightMargin(0.1);
  tdrStyle->SetPadTopMargin(0.08);
  tdrStyle->SetPadBottomMargin(0.15);
  tdrStyle->SetLegendBorderSize(0);

  const int nBins = 9;
  int bins[nBins+1] = {20,30,40,50,60,80,100,140,200,400};
  TFile* f = new TFile("FitterResults_test_Icone5.root");

  TH1F* refHist = (TH1F*)f->Get("MeanRefPt");


  TDirectory* dir = (TDirectory*)f; 
  dir->cd("FittedHistograms");

  float RefPt[nBins], errRefPt[nBins], Resp[nBins], errResp[nBins]; 


  for(int i=0; i<nBins; i++) {

    RefPt[i] = refHist->GetBinContent(i+1);
    errRefPt[i] = refHist->GetBinError(i+1);


    TString histname = Form("responseHistReco_%d_%d", bins[i], bins[i+1]);
    //  TString histname = Form("responseHistGen_%d_%d", bins[i], bins[i+1]);

    TH1F* h1 = (TH1F*)gDirectory->Get(histname);
    TF1 *myfunc = h1->GetFunction("g");

    if(myfunc==0) {
      Resp[i] = 0.0;
      errResp[i] = 0.0;
    }
    else {
      Resp[i] = myfunc->GetParameter(1);
      errResp[i] = myfunc->GetParError(1);
    }
  }


  TGraphErrors* plot = 
    new TGraphErrors( nBins, RefPt, Resp, errRefPt, errResp); 

  TAxis* xaxis = plot->GetXaxis();
  TAxis* yaxis = plot->GetYaxis();

  xaxis->SetTitle("p_{T}^{Z} (GeV/c)");
  xaxis->SetNdivisions(505);
  xaxis->SetTitleOffset(1.0);
  double xmin = xaxis->GetXmin();
  double xmax = xaxis->GetXmax();

  yaxis->SetTitle("Response");
  yaxis->SetNdivisions(505);
  yaxis->SetTitleOffset(1.4);
  plot->SetMaximum(1.04);
  plot->SetMinimum(0.96);

  TLine* line = new TLine(xmin, 1, xmax, 1);
  line->SetLineColor(2);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TCanvas* can = new TCanvas( "can", "", 500, 500);
  plot->Draw("AP");
  line->Draw();
  //  delete f;


}
