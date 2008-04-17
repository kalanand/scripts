void plotResponseMany() {

  int bins[11] = {20,30,40,50,60,80,100,140,200,400,800};
  // int bins[20] = {20, 30, 50, 80, 120, 170, 230, 300, 380, 470, 600, 800, 1000, 1400, 1800, 2200, 2600, 3000, 3500, 5000};



  // TFile* f = new TFile("MPV_FitterResults.root");
  TFile* f = new TFile("FitterResults_test_Icone5.root");
  //  TFile* f = new TFile("Corrected/FitterResults_test_Icone5.root");
  // TFile* f = new TFile("noCut/FitterResults_test_Icone5.root");

  TDirectory* dir = (TDirectory*)f; 
  dir->cd("FittedHistograms");


  for(int i=0; i<10; i++) {
    TString histname = Form("responseHistReco_%d_%d", bins[i], bins[i+1]);
    //  TString histname = Form("responseHistGen_%d_%d", bins[i], bins[i+1]);
    
        TString plotname = Form("Trunc_Uncorrected_response_%d_%d", bins[i], bins[i+1]);
    // TString plotname = Form("Trunc_Corrected_response_%d_%d", bins[i], bins[i+1]);
    // TString plotname = Form("Trunc_Generated_response_%d_%d", bins[i], bins[i+1]);
    TString label = Form("#color[4]{%d < Z p_{T} < %d GeV/c}", 
			 bins[i], bins[i+1]);
    TH1F* h1 = (TH1F*)gDirectory->Get(histname);
    makeplot(*h1, plotname, label);
  }

  delete f;

}





void makeplot(TH1& hist, TString plotname, TString label) {
  gROOT->ProcessLine(".L mystyle.C");
  setTDRStyle();
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.1);
  tdrStyle->SetPadTopMargin(0.08);
  tdrStyle->SetLegendBorderSize(0);

  TAxis* xaxis = hist.GetXaxis();
  TAxis* yaxis = hist.GetYaxis();

  xaxis->SetTitle("Response   ");
  xaxis->SetNdivisions(505);
  xaxis->SetTitleOffset(1.3);
  yaxis->SetTitle("Events   ");
  yaxis->SetNdivisions(505);
  yaxis->SetTitleOffset(1.3);


  TCanvas* can = new TCanvas(TString("can_")+hist.GetName(), "", 500, 500);
  gStyle->SetStatFormat("6.2g");
  gStyle->SetStatH(.45);
  gStyle->SetStatW(.3);
  //  gStyle->SetStatW(.2);
  gStyle->SetOptFit(1110);
  gStyle->SetOptStat(2200);


//   leg_hist = new TLegend(0.52,0.6,0.86,0.75);
//   leg_hist = new TLegend(0.3,0.92,0.6,0.998);

  leg_hist = new TLegend(0.1,0.92,0.4,0.998);
  leg_hist->AddEntry("", label,"");
  leg_hist->SetFillColor(0);
  hist.Draw("e");
  leg_hist->Draw();

  can->SaveAs( plotname+TString(".eps") );
  can->SaveAs( plotname+TString(".gif") );
  can->SaveAs( plotname+TString(".root") );
  delete can;
  delete leg_hist;

}



