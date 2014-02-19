#include "TGraphErrors.h"
#include "TFile.h"
#include "TKey.h"
#include "TObjArray.h"
#include <iostream>
#include <vector>


gROOT->Reset();
gROOT->SetStyle("Plain");
gStyle->SetOptStat(0000);
gStyle->SetOptFit(0000); 
gStyle->SetPalette(1);

const double rMin = 8.2;
const double rMax = 700.0;
//const double rMax = 2000.0;
const int rBin = (int)(rMax - rMin);




void Combiner () {
  //  Combiner ( "L3Graphs_test_Icone5.root", "csa08.root");
  Combiner ( "L3Graphs_test_Icone5.root", "summer08.root", 
	     "corrections_responses_Summer08_Zmumujet_100pb.root", 1);
}





void Combiner ( char* ZeeJetFile, char* GammaJetFile, 
		char* ZmmJetFile, const int  plotCase) {


  // Get Zee + Jet
  TFile* inf = new TFile( ZeeJetFile, "r" );
  TGraphErrors *g_Cor = (TGraphErrors*)inf->Get("Correction_vs_CaloPt");
  g_Cor->SetMarkerColor(2);
  g_Cor->SetLineColor(2);
  g_Cor->SetMarkerStyle(20);
  TMatrixD *COV_Cor = (TMatrixD*)inf->Get("CovMatrix_Correction");
  TF1 *CorFit = (TF1*)g_Cor->GetFunction("CorFit");
  CorFit->SetLineColor(2);
  CorFit->SetLineWidth(2);
  CorFit->SetRange(rMin,rMax);  
  TH1F* hCorUncertainty = Uncertainty(*CorFit, *COV_Cor, "CorrErr0");
  TF1* MCOverCor = TransformationFunction( *CorFit);
  MCOverCor->SetLineColor(2);



  // Get dijet MC curve
  TF1* MCcurve = Dijet_MC_Curve("MCcurve");



  // Get Gamma + Jet
  TFile* inf2 = new TFile( GammaJetFile, "r" );

  //////////////////////////////////////////////////////////////////
  // reverse engineering to get Mikko's "correct" correction plot !
  /////////////////////////////////////////////////////////////////

  TGraphErrors *gammaJetResp = 
    (TGraphErrors*)inf2->Get("measrespvsptphot_sig");
  TGraphErrors *gammaJetCorr = 
    (TGraphErrors*)inf2->Get("corr_sig_100pb");

  Double_t* vectX = gammaJetCorr->GetX(); 
  Double_t* vectEX = gammaJetCorr->GetEX();
  Double_t* vectY = gammaJetResp->GetY(); 
//   Double_t* vectEY = gammaJetResp->GetEY();
  Double_t* vectEY = gammaJetCorr->GetEY();

  int mumPts = gammaJetCorr->GetN();

  for (int j=0; j<mumPts; ++j) {
    vectY[j] = 1/vectY[j];
    vectEY[j] = pow(vectY[j], 2) *vectEY[j]; 
  }

  TGraphErrors *g_Cor2 = new TGraphErrors(mumPts,vectX,vectY,vectEX,vectEY);
  // TGraphErrors *g_Cor2 = (TGraphErrors*)inf2->Get("corr_sig_100pb");
  g_Cor2->SetMarkerSize(0.9);
  g_Cor2->SetMarkerColor(1);
  g_Cor2->SetLineColor(1);






  TF1* CorFit2 = new TF1("CorFit2","[0]+[1]/(pow(log10(x),[2])+[3])",
		       rMin,rMax);
  g_Cor2->Fit(CorFit2,"RQ");
  g_Cor2->SetMarkerStyle(20);
  //  TF1 *CorFit2 = (TF1*)g_Cor2->GetFunction("biassum_sig");
  //  CorFit2->SetRange(rMin,rMax);  
  TF1* MCOverCor2 = TransformationFunction( *CorFit2);
  TPaveStats *p1 
    = (TPaveStats*)g_Cor2->GetListOfFunctions()->FindObject("stats");
  g_Cor2->GetListOfFunctions()->Remove(p1);



  // Get Zmumu + Jet
  TFile* inf3 = new TFile( ZmmJetFile, "r" );
  TGraphAsymmErrors *g_Cor3 = 
    (TGraphAsymmErrors*)inf3->Get("L2CorJetIC5CaloCorrections_toys");
  g_Cor3->SetMarkerColor(7);
  g_Cor3->SetLineColor(7);
  g_Cor3->SetMarkerStyle(22);
  TF1 *CorFit3 = 
    (TF1*)g_Cor3->GetFunction("L2CorJetIC5Calo_toys_corrections_curve");
  CorFit3->SetRange(rMin,rMax);  
  TF1* MCOverCor3 = TransformationFunction( *CorFit3);
  MCOverCor3->SetLineColor(7);



  //////////////////////////////////////////////////
  //////////////////////////////////////////////////


  if(plotCase==1) 
    PlotTransFunction( hCorUncertainty, *g_Cor, *g_Cor2, *g_Cor3,
		       *MCOverCor, *MCOverCor2, *MCOverCor3);

  // Apply Transformation function
  TGraphErrors* new_gr = ScaleGraph( *g_Cor, *MCOverCor);
  TGraphErrors* new_gr2 = ScaleGraph( *g_Cor2, *MCOverCor2);
  TGraphErrors* new_gr3 = ScaleGraph( *g_Cor3, *MCOverCor3);
  new_gr->SetMarkerColor(2);
  new_gr->SetLineColor(2);

  // Now combine all the L3 corrections
  //  TGraphErrors* merged_gr = MergeGraph ( *new_gr, *new_gr2 );
  TGraphErrors* merged_gr = MergeGraph ( *new_gr, *new_gr2,  *new_gr3 );

  if(plotCase>1 && plotCase<6) 
    PlotCombinedCorrection( *merged_gr, *new_gr, 
			    *new_gr2,  *new_gr3, plotCase);
  if(plotCase==6) PlotConstantDisagreement( *merged_gr );
  if(plotCase==7) PlotLinearDisagreement( *merged_gr );
}

// ---------------------------------------------------









//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

// Plot individual corrections and Trans Functions.

void PlotTransFunction( TH1F* hCorUncertainty, 
			const TGraphErrors& g1, 
			const TGraphErrors& g2, 
			const TGraphAsymmErrors& g3,
			const TF1& f1, const TF1& f2, 
			const TF1& f3) {


  gROOT->ProcessLine(".L mystyle.C");
  setTDRStyle();

  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.08);
  tdrStyle->SetPadTopMargin(0.08);
  tdrStyle->SetLegendBorderSize(0);
  

//   tdrStyle->SetPadLeftMargin(0.2);
//   tdrStyle->SetPadRightMargin(0.1);
//   tdrStyle->SetPadTopMargin(0.08);
//   tdrStyle->SetLegendBorderSize(0);
  
  hCorUncertainty->SetMaximum(2.4);
  hCorUncertainty->SetMinimum(0.8);
  hCorUncertainty->SetMarkerColor(0);
  hCorUncertainty->SetLineColor(0);
  hCorUncertainty->SetFillColor(0);
  TAxis* xaxis = hCorUncertainty->GetXaxis();
  TAxis* yaxis = hCorUncertainty->GetYaxis();
  xaxis->SetNdivisions(508);  
  xaxis->SetTitleOffset(1.34);
  xaxis->SetRangeUser(5,800);
  xaxis->SetMoreLogLabels();
  xaxis->SetNoExponent();
  yaxis->SetNdivisions(414);  
  yaxis->SetTitleOffset(1.4);
  yaxis->SetRangeUser(0.8,3.0);


  // if we don't want to draw the fit curves
  TF1 *cc1 = (TF1*)g1.GetFunction("CorFit");
  TF1 *cc2 = (TF1*)g2.GetFunction("CorFit2");
  TF1 *cc3 = (TF1*)g3.GetFunction("L2CorJetIC5Calo_toys_corrections_curve");
  cc1->SetParameter(0, 100);
  cc2->SetParameter(0, 100);
  cc3->SetParameter(0, 100);
    

  // Get dijet MC curve
  TF1* MCcurve = Dijet_MC_Curve("MCcurve");
  MCcurve->SetLineWidth(1);

  TCanvas *c_Correction = new TCanvas("Correction","Correction",500,500);
  c_Correction->cd(); 
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gPad->SetLogx();
  g1.RemovePoint(9);
  g3.RemovePoint(0);

  hCorUncertainty->Draw("E3");
  g1.Draw("Psame"); 
  g2.Draw("Psame"); 
  g3.Draw("Psame"); 
  f1.Draw("Lsame"); 
  f2.SetLineColor(1);
  f2.Draw("Lsame"); 
  f3.Draw("Lsame"); 
  MCcurve->Draw("Lsame"); 

  TLegend *leg = new TLegend(0.45,0.5,0.89,0.88);
  leg->AddEntry( &g1,"Correction: Z#rightarrow e^{+}e^{-}+jet","P");
  leg->AddEntry( &g2,"Correction: #gamma+jet","P");
  leg->AddEntry( &g3,"Correction: Z#rightarrow #mu^{+}#mu^{-}+jet","le");
  leg->AddEntry(MCcurve,"Correction: Dijet MC-truth","L");

  leg->AddEntry( &f1,"Trans. function: Zee+jet","L");
  leg->AddEntry( &f2,"Trans. function: #gamma+jet","L");
  leg->AddEntry( &f3,"Trans. function: Z#mu#mu+jet","L");

  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetMargin(0.15);
  leg->Draw();
  c_Correction->SaveAs("CombinedCorrection-1.gif");
  c_Correction->SaveAs("CombinedCorrection-1.eps");
  c_Correction->SaveAs("CombinedCorrection-1.root");
}




//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

// Plot combined correction
// "option" can be 2,3,4,5

void PlotCombinedCorrection( const TGraphErrors& g,  
			     const TGraphErrors& g1, 
			     const TGraphErrors& g2, 
			     const TGraphErrors& g3, const int option) {

  gROOT->ProcessLine(".L mystyle.C");
  setTDRStyle();
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.08);
  tdrStyle->SetPadTopMargin(0.08);
  tdrStyle->SetLegendBorderSize(0);
  
  TF1* C_all = new TF1("C_all","[0]+[1]/(pow(log10(x),[2])+[3])",rMin,rMax);
  g.Fit(C_all,"RQ");
  fitter = TVirtualFitter::GetFitter();
  TMatrixD* cov = new TMatrixD(4,4,fitter->GetCovarianceMatrix());
  TH1F* hCorrErr = Uncertainty(*C_all, *cov, "CorrErr");

  C_all->SetLineColor(6);
  C_all->SetMarkerColor(6);
  C_all->SetLineWidth(2);
  if(option==4) hCorrErr->SetMaximum(4);
  else hCorrErr->SetMaximum(2.8);
  TAxis* xaxis = hCorrErr->GetXaxis();
  xaxis->SetNdivisions(508);
  xaxis->SetNdivisions(508);
  if(option==4) xaxis->SetRangeUser(1,2000);
  else xaxis->SetRangeUser(8,800);
  xaxis->SetMoreLogLabels();
  xaxis->SetNoExponent();
  hCorrErr->GetYaxis()->SetTitleOffset(1.2);
  if(option==4) xaxis->SetTitleOffset(1.2);
  else xaxis->SetTitleOffset(1);
  g1.SetMarkerStyle(20);
  g2.SetMarkerStyle(20);
  g2.SetMarkerColor(1);
  g1.SetMarkerSize(1.2);
  g2.SetMarkerSize(0.9);
  g2.SetLineColor(1);
  g3.SetMarkerStyle(20);
  g3.SetMarkerColor(7);
  g3.SetMarkerSize(0.9);
  g3.SetLineColor(7);

  TF1* MCcurve = Dijet_MC_Curve("MCcurve");
  tdrStyle->SetPadTopMargin(0.08);
  tdrStyle->SetPadBottomMargin(0.16);

  if(option==4) { 
    MCcurve->SetRange(1.2, 2000);
    C_all->SetRange(1.2, 2000);
  }
  g1.RemovePoint(9);
  g3.RemovePoint(0);

  TCanvas *c_Correction3 = new TCanvas("Correction3","Correction3",500,500);
  gStyle->SetOptStat(0000);
  c_Correction3->cd(); 
  hCorrErr->Draw("E3");
  g3.Draw("P");
  g2.Draw("P"); 
  g1.Draw("P"); 
  if (option>3) MCcurve->SetLineWidth(2);
  if (option>2) MCcurve->Draw("Lsame");
  C_all->Draw("Lsame");
  TLegend *leg = new TLegend(0.4,0.55,0.85,0.85);
  leg->AddEntry( &g1,"Z#rightarrow e^{+}e^{-}+jet","LP");
  leg->AddEntry( &g3,"Z#rightarrow #mu^{+}#mu^{-}+jet","LP");
  leg->AddEntry( &g2,"#gamma+jet","LP");
  leg->AddEntry(C_all,"Combined correction","LP");
  leg->AddEntry(hCorrErr,"Combined uncertainty","F");
  if (option>2) leg->AddEntry(MCcurve,"Dijet MC-truth","L");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
  gPad->SetLogx();
  //  if(option==4) 
  TString plotname = Form("CombinedCorrection-%d", option);
  c_Correction3->SaveAs( plotname + TString(".gif") );
  c_Correction3->SaveAs( plotname + TString(".eps") );
  c_Correction3->SaveAs( plotname + TString(".root") );
}



//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

 // Assume that the combined L3 corrections points in data are 10% higher 

void PlotConstantDisagreement( const TGraphErrors& g ) {

  TGraphErrors* gr_10pc = ScaleGraph( g, 1.1);
  TF1* C_10pc = new TF1("C_10pc","[0]+[1]/(pow(log10(x),[2])+[3])",rMin,rMax);
  gr_10pc->Fit(C_10pc,"RQ");
  fitter = TVirtualFitter::GetFitter();
  TMatrixD* cov_10pc = new TMatrixD(4,4,fitter->GetCovarianceMatrix());
  TH1F* hCorrErr_10pc = Uncertainty(*C_10pc, *cov_10pc, "CorrErrAll_10pc");
  C_10pc->SetLineColor(6);

  TF1* MCcurve = Dijet_MC_Curve();

  TCanvas *c_Correction4 = new TCanvas("Correction4","Correction4",500,500);
  gStyle->SetOptStat(0000);
  c_Correction4->cd(); 
  gPad->SetLogx();  
  hCorrErr_10pc->SetMaximum(4.8);
  hCorrErr_10pc->Draw("E3");
  gr_10pc->Draw("Psame"); 
  C_10pc->Draw("Lsame");
  MCcurve->Draw("Lsame");
  TLegend *leg = new TLegend(0.5,0.65,0.89,0.89);
  leg->AddEntry(gr_10pc,"Measurement from data","LP");
  leg->AddEntry(C_10pc,"Curve from fitting the data points","LP");
  leg->AddEntry(MCcurve,"Curve from dijet MC-truth","L");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
}





//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////


// Assume that the combined L3 corrections points in data are 10% higher 
// at lowest Pt and only 2% higher at highest pT

void PlotLinearDisagreement( const TGraphErrors& g ) {

  gROOT->ProcessLine(".L mystyle.C");
  setTDRStyle();
  tdrStyle->SetPadLeftMargin(0.2);
  tdrStyle->SetPadRightMargin(0.1);
  tdrStyle->SetPadTopMargin(0.08);
  tdrStyle->SetLegendBorderSize(0);


  TF1* linear = new TF1("linear","1.1613 - 0.0471 * log10(x)",rMin,rMax);
  TGraphErrors* gr_linear = ScaleGraph( g, *linear);
  TF1* C_linear = new TF1("C_linear","[0]+[1]/(pow(log10(x),[2])+[3])",
			  rMin,rMax);
  gr_linear->Fit(C_linear,"RQ");
  fitter = TVirtualFitter::GetFitter();
  TMatrixD* cov_linear = new TMatrixD(4,4,fitter->GetCovarianceMatrix());
  TH1F* hCorrErr_linear = Uncertainty(*C_linear, *cov_linear, "CorrErr");
  C_linear->SetLineColor(6);

  TF1* MCcurve = Dijet_MC_Curve();


  TCanvas *c_Correction5 = new TCanvas("Correction5","Correction5",500,500);
  gStyle->SetOptStat(0000);
  c_Correction5->cd(); 
  gPad->SetLogx();  
  hCorrErr_linear->SetMaximum(5.05);
  hCorrErr_linear->Draw("E3");
  gr_linear->Draw("Psame"); 
  C_linear->Draw("Lsame");
  MCcurve->Draw("Lsame");
  TLegend *leg = new TLegend(0.5,0.65,0.89,0.89);
  leg->AddEntry(gr_linear,"Measurement from data","LP");
  leg->AddEntry(C_linear,"Curve from fitting the data points","LP");
  leg->AddEntry(MCcurve,"Curve from dijet MC-truth","L");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
}



//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

double FitUncertainty( TF1* f, TMatrixD* COV, double x) {
  int i,j;
  double df,sum,y,z,x;
  double PartialDerivative[10],Parameter[10];
  int npar = 4;
  int N = f->GetNumberFreeParameters();
  int dim = COV->GetNrows();

  if (dim != npar || N != npar)
    {
      cout<<"ERROR: wrong number of parameters !!!!"<<endl;
      return(-1);
    }  
  for(i=0;i<npar;i++)
    Parameter[i] = f->GetParameter(i);
  z = pow(log10(x),Parameter[2]);  
  PartialDerivative[0] = 1.;
  PartialDerivative[1] = 1./(z+Parameter[3]);
  PartialDerivative[3] = -Parameter[1]/pow(z+Parameter[3],2);
  PartialDerivative[2] = PartialDerivative[3]*log(log10(x))*z;
  sum = 0.;
  for(i=0;i<npar;i++)
    for(j=0;j<npar;j++)
      {
        y = PartialDerivative[i]*PartialDerivative[j]*COV(i,j);
        sum+=y;
      }
  df = sqrt(sum);
  return df;
}





//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

TGraphErrors* MergeGraph(const TGraphErrors& g1, 
			 const TGraphErrors& g2, 
			 const TGraphErrors& g3) {

  TGraphErrors* temp = MergeGraph ( g1, g2 );
  TGraphErrors* g = MergeGraph ( *temp, g3 );
  delete temp;
  return g;
}


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

TGraphErrors* MergeGraph(const TGraphErrors& g1, 
			 const TGraphErrors& g2) {
  vector<double> x;
  vector<double> y;
  vector<double> xerr;
  vector<double> yerr;
 
  for ( int i = 0; i < g1.GetN(); ++i ) {
    x.push_back(g1.GetX()[i]);
    xerr.push_back(g1.GetEX()[i]);
    y.push_back(g1.GetY()[i]);
    yerr.push_back(g1.GetEY()[i]);
  }

  for ( int i = 0; i < g2.GetN(); ++i ) {
    x.push_back(g2.GetX()[i]);
    xerr.push_back(g2.GetEX()[i]);
    y.push_back(g2.GetY()[i]);
    yerr.push_back(g2.GetEY()[i]);
  }
  
  TGraphErrors* g = new TGraphErrors(x.size(),
				     &(x[0]),&(y[0]),
                                     &(xerr[0]),&(yerr[0]));
  g->SetName(g1.GetName());
  return g;  
}




//////////////////////////////////////////////////
//////////////////////////////////////////////////
TGraphErrors* ScaleGraph(const TGraphAsymmErrors& g, 
			 const TF1& f) {

  vector<double> x;
  vector<double> y;
  vector<double> xerr;
  vector<double> yerr;

  for(int i=0; i<g.GetN(); ++i) {
    double xi = g.GetX()[i];
    x.push_back(xi);
    xerr.push_back( g.GetEXhigh()[i] );
    y.push_back( g.GetY()[i] * f.Eval(xi) );
    yerr.push_back( g.GetEYhigh()[i] * f.Eval(xi) );
  }

  TGraphErrors* new_gr = new TGraphErrors(x.size(),
				     &(x[0]),&(y[0]),
                                     &(xerr[0]),&(yerr[0]));
  new_gr->SetName(g.GetName());
  new_gr->SetTitle("");
  new_gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  new_gr->GetYaxis()->SetTitle("Correction factor");
  return new_gr;  
}




//////////////////////////////////////////////////
//////////////////////////////////////////////////
TGraphErrors* ScaleGraph(const TGraphErrors& g, 
			 const TF1& f) {

  vector<double> x;
  vector<double> y;
  vector<double> xerr;
  vector<double> yerr;

  for(int i=0; i<g.GetN(); ++i) {
    double xi = g.GetX()[i];
    x.push_back(xi);
    xerr.push_back( g.GetEX()[i] );
    y.push_back( g.GetY()[i] * f.Eval(xi) );
    yerr.push_back( g.GetEY()[i] * f.Eval(xi) ); 
  }

  TGraphErrors* new_gr = new TGraphErrors(x.size(),
				     &(x[0]),&(y[0]),
                                     &(xerr[0]),&(yerr[0]));
  new_gr->SetName(g.GetName());
  new_gr->SetTitle("");
  new_gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  new_gr->GetYaxis()->SetTitle("Correction factor");
  return new_gr;  
}



//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
TGraphErrors* ScaleGraph(const TGraphErrors& g, 
			 const double f) {

  vector<double> x;
  vector<double> y;
  vector<double> xerr;
  vector<double> yerr;

  for(int i=0; i<g.GetN(); ++i) {
    double xi = g.GetX()[i];
    x.push_back(xi);
    xerr.push_back( g.GetEX()[i] );
    y.push_back( g.GetY()[i] * f );
    yerr.push_back( g.GetEY()[i] * f ); 
  }

  TGraphErrors* new_gr = new TGraphErrors(x.size(),
				     &(x[0]),&(y[0]),
                                     &(xerr[0]),&(yerr[0]));
  new_gr->SetName(g.GetName());
  new_gr->SetTitle("");
  new_gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  new_gr->GetYaxis()->SetTitle("Correction factor");
  return new_gr;  
}


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
TF1* TransformationFunction(const TF1& f) {
  
  // Get dijet MC curve
  TF1* f0 = Dijet_MC_Curve("f0");

  // This is from X+jet MC
  TF1 f1("f1","[0]+[1]/(pow(log10(x),[2])+[3])", rMin, rMax); 
  f1.SetParameter(0, f.GetParameter(0));
  f1.SetParameter(1, f.GetParameter(1));
  f1.SetParameter(2, f.GetParameter(2));
  f1.SetParameter(3, f.GetParameter(3));
  
  // take ratio
  TF1* tf = new TF1("tf","f0/f1", rMin, rMax); 
  tf->SetLineColor(2);
  tf->SetLineWidth(2);

  tf->SetName(f.GetName());
  return tf;  
}


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

TH1F* Uncertainty(const TF1& f, TMatrixD& COV, 
		  char* name, Color_t color=5) {
  
  TH1F* h = new TH1F( name, "", 100010, 1, 20001);
  double x,y,e;

  for(int i=0; i<100010; i++) {
      x = h->GetBinCenter(i+1);
      y = f.Eval(x);
      e = FitUncertainty( &f, &COV, x );
      h->SetBinContent( i+1, y);
      h->SetBinError( i+1, e ); 
    }

  h->SetMinimum(1);
  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h->GetYaxis()->SetTitle("Correction factor");  
  h->SetLineColor(5);
  h->SetFillColor(color);
  h->SetMarkerColor(color);  

  return h;  
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

TF1* Dijet_MC_Curve(TString name) {
  
  // TF1* C = new TF1("C","[0]+[1]/(pow(log10(x),[2])+[3])", rMin, rMax); 
TF1* C = new TF1("C","[0]+[1]/(pow(log10(x),[2])+[3])", 1.2, rMax); 

  /////// for CSA07 correction
  //   C->SetParameter(0, 0.996998);
  //   C->SetParameter(1, 4.39412);
  //   C->SetParameter(2, 2.96134);
  //   C->SetParameter(3, 1.69966);

  /////// for Summer08 correction
  C->SetParameter(0, 0.998293);
  C->SetParameter(1, 5.43056);
  C->SetParameter(2, 3.3444);
  C->SetParameter(3, 2.39809);


  C->SetLineColor(4);
  // C->SetLineStyle(2);
  C->SetLineWidth(4);
  C->SetName(name);
  return C;
} 

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
