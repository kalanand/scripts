{
gROOT->ProcessLine(".L mystyle.C");
setTDRStyle();
tdrStyle->SetPadLeftMargin(0.16);
tdrStyle->SetPadRightMargin(0.1);
tdrStyle->SetPadTopMargin(0.08);
tdrStyle->SetLegendBorderSize(0);


const int NGenPtBins=21;
const double GenPt[NGenPtBins+1] = {0,15,20,30,50,80,120,170,230,300,380,470,600,800,1000,1400,1800,2200,2600,3000,3500,-1};
TString ptBin[NGenPtBins];
TString basename = "/uscms_data/d1/kalanand/trash/CSA07-ZeeJets_Pt_";
TFile* f[NGenPtBins];
const double lumi = 100.0;



double crosssection[NGenPtBins];
crosssection[0]  =  104400.0 * 0.98;
crosssection[1]  =  353.1          ;
crosssection[2]  =  326.7          ;
crosssection[3]  =  227.0          ;
crosssection[4]  =  93.17          ;
crosssection[5]  =  31.48          ;
crosssection[6]  =  9.63           ;
crosssection[7]  =  2.92           ;
crosssection[8]  =  0.8852         ;
crosssection[9]  =  0.2936         ;
crosssection[10] =  0.1025         ;
crosssection[11] =  0.04242        ;
crosssection[12] =  0.01443        ;
crosssection[13] =  0.002859       ;
crosssection[14] =  0.00094        ;
crosssection[15] =  0.00009536     ;
crosssection[16] =  0.00001232     ;
crosssection[17] =  0.000001839    ;
crosssection[18] =  0.0000002881   ;
crosssection[19] =  0.00000004764  ;
crosssection[20] =  0.000000004515 ;


for(int i=0; i<NGenPtBins; i++) {
  ptBin[i] = Form("%d_%d", (int) GenPt[i], (int) GenPt[i+1] );
  f[i] = new TFile( basename+ptBin[i]+TString(".root") );
}



TProfile normResidual_Zpt("normResidual_Zpt", "", 20, 0, 500, -0.12, 0.12);
TProfile normResZptGen("normResZptGen","",20, 0, 500, -0.1, 0.1);
TProfile normResidual_empt("normResidual_empt", "", 20, 0, 200, -0.12, 0.12);
TProfile normResidual_emptGen("normResidual_emptGen", "", 20, 0, 200, -0.12, 0.12);
TProfile normResidual_ZptCut("normResidual_ZptCut", "", 20, 0, 500, -0.12, 0.12);
TProfile normResZptGenCut("normResZptGenCut","",20, 0, 500, -0.1, 0.1);
TProfile normResidual_emptCut("normResidual_emptCut", "", 20, 0, 200, -0.12, 0.12);
TProfile normResidual_emptGenCut("normResidual_emptGenCut", "", 20, 0, 200, -0.12, 0.12);



normResidual_ZptCut.SetLineColor(2);
normResZptGenCut.SetLineColor(2);
normResidual_emptCut.SetLineColor(2);
normResidual_emptGenCut.SetLineColor(2);

normResidual_ZptCut.SetMarkerColor(2);
normResZptGenCut.SetMarkerColor(2);
normResidual_emptCut.SetMarkerColor(2);
normResidual_emptGenCut.SetMarkerColor(2);



TAxis* xnr_Zpt = normResidual_Zpt.GetXaxis();
TAxis* ynr_Zpt = normResidual_Zpt.GetYaxis();
xnr_Zpt->SetTitle("Reconstructed Z p_{T} (GeV/c)");
xnr_Zpt->SetTitleSize(0.04);
xnr_Zpt->SetTitleOffset(1.4);
xnr_Zpt->SetNdivisions(505);
ynr_Zpt->SetTitle("#frac{Generated Z p_{T} - Reconstructed Z p_{T}}{Generated Z p_{T}}   ");
ynr_Zpt->SetTitleSize(0.04);
ynr_Zpt->SetTitleOffset(2.4);


TAxis* xnr_ZptGen = normResZptGen.GetXaxis();
TAxis* ynr_ZptGen = normResZptGen.GetYaxis();
xnr_ZptGen->SetTitle("Generated Z p_{T} (GeV/c)");
xnr_ZptGen->SetTitleSize(0.04);
xnr_ZptGen->SetTitleOffset(1.4);
xnr_ZptGen->SetNdivisions(505);
ynr_ZptGen->SetTitle("#frac{Generated Z p_{T} - Reconstructed Z p_{T}}{Generated Z p_{T}}   ");
ynr_ZptGen->SetTitleSize(0.04);
ynr_ZptGen->SetTitleOffset(2.4);


TAxis* xnr_empt = normResidual_empt.GetXaxis();
TAxis* ynr_empt = normResidual_empt.GetYaxis();
xnr_empt->SetTitle("Reconstructed e^{-} p_{T} (GeV/c)");
xnr_empt->SetTitleSize(0.04);
xnr_empt->SetTitleOffset(1.4);
xnr_empt->SetNdivisions(505);
ynr_empt->SetTitleSize(0.04);
ynr_empt->SetTitleOffset(2.5);
ynr_empt->SetTitle("#frac{Generated e^{-} p_{T} - Reconstructed e^{-} p_{T}}{Generated e^{-} p_{T}}   ");


TAxis* xnr_emptGen = normResidual_emptGen.GetXaxis();
TAxis* ynr_emptGen = normResidual_emptGen.GetYaxis();
xnr_emptGen->SetTitle("Generated e^{-} p_{T} (GeV/c)");
xnr_emptGen->SetTitleSize(0.04);
xnr_emptGen->SetTitleOffset(1.4);
xnr_emptGen->SetNdivisions(505);
ynr_emptGen->SetTitleSize(0.04);
ynr_emptGen->SetTitleOffset(2.5);
ynr_emptGen->SetTitle("#frac{Generated e^{-} p_{T} - Reconstructed e^{-} p_{T}}{Generated e^{-} p_{T}}   ");






for(int i=0; i<NGenPtBins; i++) {
  TTree* t = (TTree*) f[i].Get("ZJet");
  int nEvents  = (int) t.GetEntries();
  double coeff = lumi * crosssection[i]  /  nEvents;
  TString wt = Form("%f", coeff);  
  if(i==1) wt = Form("%f*(JetGenPt[0][0]<500.0)", coeff); 


  TString cutZCut =  wt+"*(abs(mZee-91.2)<10.0 && abs(Z_Pt/Z_PtGen)<100 && ePlusPt>20.0 && eMinusPt>20.0 && ((fabs(ePlusEta)<1.4442) || (fabs(ePlusEta)>1.560 && fabs(ePlusEta)<2.5)) && ((fabs(eMinusEta)<1.4442) || (fabs(eMinusEta)>1.560 && fabs(eMinusEta)<2.5)))";

  TString cutZ =  wt+"*(abs(mZee-91.2)<10.0 && abs(Z_Pt/Z_PtGen)<100 && ((fabs(ePlusEta)<1.4442) || (fabs(ePlusEta)>1.560 && fabs(ePlusEta)<2.5)) && ((fabs(eMinusEta)<1.4442) || (fabs(eMinusEta)>1.560 && fabs(eMinusEta)<2.5)))";


  t->Draw("(Z_PtGen-Z_Pt)/Z_PtGen:Z_Pt>>+normResidual_Zpt", cutZ, "goff");
  t->Draw("(Z_PtGen-Z_Pt)/Z_PtGen:Z_PtGen>>+normResZptGen", cutZ, "goff");
  t->Draw("(eMinusPtGen-eMinusPt)/eMinusPtGen:eMinusPt>>+normResidual_empt", cutZ, "goff");
  t->Draw("(eMinusPtGen-eMinusPt)/eMinusPtGen:eMinusPtGen>>+normResidual_emptGen", cutZ, "goff");



  t->Draw("(Z_PtGen-Z_Pt)/Z_PtGen:Z_Pt>>+normResidual_ZptCut", cutZCut, "goff");
  t->Draw("(Z_PtGen-Z_Pt)/Z_PtGen:Z_PtGen>>+normResZptGenCut", cutZCut, "goff");
  t->Draw("(eMinusPtGen-eMinusPt)/eMinusPtGen:eMinusPt>>+normResidual_emptCut", cutZCut, "goff");
  t->Draw("(eMinusPtGen-eMinusPt)/eMinusPtGen:eMinusPtGen>>+normResidual_emptGenCut", cutZCut, "goff");

}





///////////////////////////////////////
/////////////////////////////////////



TCanvas can7("can7", "", 600, 600);
gStyle->SetOptStat(0);
normResidual_Zpt->Draw("e");
normResidual_ZptCut->Draw("esame");
leg_hist = new TLegend(0.5,0.75,0.9,0.85);
leg_hist->AddEntry(&normResidual_Zpt,"No cut on electron p_{T}","l");
leg_hist->AddEntry(&normResidual_ZptCut,"electron p_{T} > 20 GeV/c","l");
leg_hist->SetFillColor(0);
leg_hist->Draw();
TString plotname = "NormRes_Zpt_20GeVcut";
can7.SaveAs( plotname+TString(".eps") );
can7.SaveAs( plotname+TString(".gif") );
can7.SaveAs( plotname+TString(".root") );



TCanvas can13("can13", "", 600, 600);
gStyle->SetOptStat(0);
normResidual_empt->Draw("e");
normResidual_emptCut->Draw("esame");
leg_hist = new TLegend(0.5,0.75,0.9,0.85);
leg_hist->AddEntry(&normResidual_empt,"No cut on electron p_{T}","l");
leg_hist->AddEntry(&normResidual_emptCut,"electron p_{T} > 20 GeV/c","l");
leg_hist->SetFillColor(0);
leg_hist->Draw();
TString plotname = "NormRes_empt_20GeVcut";
can13.SaveAs( plotname+TString(".eps") );
can13.SaveAs( plotname+TString(".gif") );
can13.SaveAs( plotname+TString(".root") );



TCanvas can19("can19", "", 600, 600);
gStyle->SetOptStat(0);
normResZptGen->Draw("e");
normResZptGenCut->Draw("esame");
leg_hist = new TLegend(0.5,0.75,0.9,0.85);
leg_hist->AddEntry(&normResZptGen,"No cut on electron p_{T}","l");
leg_hist->AddEntry(&normResZptGenCut,"electron p_{T} > 20 GeV/c","l");
leg_hist->SetFillColor(0);
leg_hist->Draw();
TString plotname = "NormRes_ZptGen_20GeVcut";
can19.SaveAs( plotname+TString(".eps") );
can19.SaveAs( plotname+TString(".gif") );
can19.SaveAs( plotname+TString(".root") );



TCanvas can23("can23", "", 600, 600);
gStyle->SetOptStat(0);
normResidual_emptGen->Draw("e");
normResidual_emptGenCut->Draw("esame");
leg_hist = new TLegend(0.5,0.75,0.9,0.85);
leg_hist->AddEntry(&normResidual_emptGen,"No cut on electron p_{T}","l");
leg_hist->AddEntry(&normResidual_emptGenCut,"electron p_{T} > 20 GeV/c","l");
leg_hist->SetFillColor(0);
leg_hist->Draw();
TString plotname = "NormRes_emptGen_20GeVcut";
can23.SaveAs( plotname+TString(".eps") );
can23.SaveAs( plotname+TString(".gif") );
can23.SaveAs( plotname+TString(".root") );

}
