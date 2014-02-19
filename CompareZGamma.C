#include "TGraphErrors.h"
#include "TGraph.h"
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

const double rMin = 20.0;
const double rMax = 820.0;
const int rBin = (int)(rMax - rMin);


void CompareZGamma() {
const int NPtBins=15;
const int NETA = 1;
 const double Pt[NPtBins+1] = {25,30,36,43,51,61,73,87,104,124,148,
				   177,212,254,304,364};
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
  TGraphErrors *respZee = new TGraphErrors(N,xRefPt,
						 yRefPt,exRefPt,eyRefPt);
  respZee->SetMarkerStyle(20);
  respZee->SetMarkerSize(1.7);
  respZee->SetLineColor(2);
  respZee->SetMarkerColor(2);

  gROOT->ProcessLine(".L mystyle.C");
  setTDRStyle();
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPadLeftMargin(0.2);
  tdrStyle->SetPadRightMargin(0.10);
  tdrStyle->SetPadBottomMargin(0.16);
  tdrStyle->SetLegendBorderSize(0);
  tdrStyle->SetTitleYOffset(1.3);



  //TFile fZmm("corrections_responses_Summer08_Zmumujet_100pb.root","read");
//   TFile* fZmm = new TFile("responses_Summer08_Zmumujet_MCerrors/Response_MC_errors.root");  
//   TGraphAsymmErrors *rZmm = 
//     (TGraphAsymmErrors *) fZmm->Get("iterativeCone5CaloJets"); 


  TGraphAsymmErrors *grae = new TGraphAsymmErrors(15);
  grae->SetName("hm_iterativeCone5CaloJets_ratio_graph");
  grae->SetTitle("hm_iterativeCone5CaloJets <Pt_{Jet}/Pt_{Z}> vs Pt_{Z}");
  grae->SetLineWidth(2);
  grae->SetMarkerStyle(3);
  grae->SetPoint(0,27.66014,0.402234);
  grae->SetPointError(0,0,0,0.007289587,0.007289587);
  grae->SetPoint(1,33.06504,0.4601342);
  grae->SetPointError(1,0,0,0.005707689,0.005707689);
  grae->SetPoint(2,39.45714,0.4832026);
  grae->SetPointError(2,0,0,0.005461802,0.005461802);
  grae->SetPoint(3,46.91662,0.5039087);
  grae->SetPointError(3,0,0,0.00505666,0.00505666);
  grae->SetPoint(4,55.79232,0.5504526);
  grae->SetPointError(4,0,0,0.004945873,0.004945873);
  grae->SetPoint(5,66.6575,0.5866128);
  grae->SetPointError(5,0,0,0.004700836,0.004700836);
  grae->SetPoint(6,79.62225,0.6163132);
  grae->SetPointError(6,0,0,0.005049584,0.005049584);
  grae->SetPoint(7,94.87024,0.6479203);
  grae->SetPointError(7,0,0,0.005435905,0.005435905);
  grae->SetPoint(8,113.0874,0.6723784);
  grae->SetPointError(8,0,0,0.006032537,0.006032537);
  grae->SetPoint(9,134.8582,0.7016644);
  grae->SetPointError(9,0,0,0.007120537,0.007120537);
  grae->SetPoint(10,160.9029,0.7289301);
  grae->SetPointError(10,0,0,0.008098604,0.008098604);
  grae->SetPoint(11,192.5925,0.7430468);
  grae->SetPointError(11,0,0,0.01016994,0.01016994);
  grae->SetPoint(12,230.49,0.7636903);
  grae->SetPointError(12,0,0,0.01296528,0.01296528);
  grae->SetPoint(13,275.7171,0.7798821);
  grae->SetPointError(13,0,0,0.01646034,0.01646034);
  grae->SetPoint(14,329.9,0.7930674);
  grae->SetPointError(14,0,0,0.02150947,0.02150947);


  TGraphErrors *respZmm = ConvertFromAsymGraph(*grae);
  respZmm->SetMarkerColor(kGreen+2);
  respZmm->SetLineColor(kGreen+2);
  respZmm->SetMarkerSize(2.0);
  respZmm->SetMarkerStyle(22);


  TFile* fGamma = new TFile("summer08.root","read"); 
  TGraphErrors *respGamma = 
    (TGraphErrors*)fGamma->Get("measrespvsptphot_sig");
  respGamma->SetMarkerColor(1);
  respGamma->SetLineColor(1);
  respGamma->SetMarkerStyle(24);
  respGamma->SetMarkerSize(1.6);
  TF1 *cc1 = (TF1*)respGamma->GetFunction("fsig");
  delete cc1;
  TPaveStats *p1 
    = (TPaveStats*) respGamma->GetListOfFunctions()->FindObject("stats");
  respGamma->GetListOfFunctions()->Remove(p1);



  // Now combine all the responses into a single response 
  TGraphErrors* gr = MergeGraph ( *respZee, *respZmm,  *respGamma );
  

  // Fit the combined response 
  TF1* fit = new TF1("fit","[0]-[1]/(pow(log10(x),[2])+[3])+[4]/x",
		       rMin,rMax);
  fit->SetParameter(0,1.);
  fit->SetParameter(1,1.);
  fit->SetParameter(2,1.);
  fit->SetParameter(3,1.);
  fit->SetParameter(4,1.);
  fit->SetLineColor(6);
  fit->SetLineWidth(2);
  gr->Fit(fit,"RQ");
  fitter = TVirtualFitter::GetFitter();
  TMatrixD* cov = new TMatrixD(5,5,fitter->GetCovarianceMatrix());
  TH1F* hCorrErr = Uncertainty(*fit, *cov, "CorrErr");
  
  TAxis* xaxis = hCorrErr->GetXaxis();
  TAxis* yaxis = hCorrErr->GetYaxis();
  xaxis->SetTitle("p_{T}^{Z/#gamma} (GeV/c)");
  yaxis->SetTitle("p_{T}^{jet} / p_{T}^{Z/#gamma}   ");
  xaxis->SetMoreLogLabels();
  xaxis->SetNoExponent();
  xaxis->SetLabelSize(0.06);
  xaxis->SetTitleSize(0.06);
  xaxis->SetTitleOffset(1.1);
  yaxis->SetLabelSize(0.06);
  yaxis->SetTitleSize(0.06);
  yaxis->SetTitleOffset(1.5);
  yaxis->SetNdivisions(510);
  hCorrErr->SetMaximum(1);
  hCorrErr->SetMinimum(0.3);


  TCanvas* Response = new TCanvas("Response","Response",500,500);
  gStyle->SetOptStat(0000);
  hCorrErr->Draw("E3");
  respGamma->Draw("P");
  respZee->Draw("P");
  respZmm->Draw("P");
  fit->Draw("Lsame");
  TLegend* leg = new TLegend(0.23,0.7,0.63,0.9);
  leg->AddEntry( respZee,"Z (#rightarrow e^{+}e^{-}) + jet","P");
  leg->AddEntry( respZmm,"Z (#rightarrow#mu^{+}#mu^{-}) + jet","P");
  leg->AddEntry( respGamma,"#gamma + jet","P");
  leg->AddEntry( fit,"Combined response","L");
  leg->AddEntry( hCorrErr,"Combined uncertainty","F");
  leg->SetMargin(0.2);
  leg->SetFillColor(0);
  leg->Draw();
  TLatex* CMS = new  TLatex();
  CMS->SetTextAlign(12);
  CMS->SetTextSize(0.04);
  CMS->SetNDC();
  CMS->DrawLatex(0.5, 0.45, "CMS Preliminary");
  // CMS->DrawLatex(0.5, 0.35, "#intLdt = 100 pb^{-1}");
  CMS->DrawLatex(0.52, 0.4, "MC Statistics");
  gPad->SetLogx();
  Response->SaveAs("Fig6a.eps");
  Response->SaveAs("Fig6a.gif");
  Response->SaveAs("Fig6a.root");
  //   Response.Close();
  //   delete leg;
  //   delete respZee;


  Float_t x0[2] = {20, 820};
  Float_t y0[2] = {0, 0};
  TGraph* axisGr = new TGraph( 2, x0, y0 );
  TAxis* xax = axisGr->GetXaxis();
  TAxis* yax = axisGr->GetYaxis();
  xax->SetTitle("p_{T}^{Z/#gamma} (GeV/c)");
  xax->SetMoreLogLabels();
  xax->SetNoExponent();
  xax->SetLabelSize(0.06);
  xax->SetTitleSize(0.06);
  xax->SetTitleOffset(1.1);
  yax->SetLabelSize(0.06);
  yax->SetTitleSize(0.06);
  yax->SetTitleOffset(1.95);
  yax->SetNdivisions(505);
  yax->SetTitle("#frac{Response - Combined Response}{Combined Response}");
  axisGr->SetMaximum(0.1);
  axisGr->SetMinimum(-0.1);
  axisGr->SetLineColor(0);
  axisGr->SetMarkerColor(0);


  TGraphErrors* nrZee = GetNormResidual( *respZee, *fit);
  nrZee->SetMarkerStyle(20);
  nrZee->SetMarkerSize(1.7);
  nrZee->SetLineColor(2);
  nrZee->SetMarkerColor(2);

  TGraphErrors* nrZmm = GetNormResidual( *respZmm, *fit);
  nrZmm->SetMarkerColor(kGreen+2);
  nrZmm->SetLineColor(kGreen+2);
  nrZmm->SetMarkerSize(2.0);
  nrZmm->SetMarkerStyle(22);

  TGraphErrors* nrGamma = GetNormResidual( *respGamma, *fit);
  nrGamma->SetMarkerColor(1);
  nrGamma->SetLineColor(1);
  nrGamma->SetMarkerStyle(24);
  nrGamma->SetMarkerSize(1.6);

  tdrStyle->SetPadLeftMargin(0.26);

  TCanvas* NormRes = new TCanvas("NormRes","NormRes",500,500);
  gStyle->SetOptStat(0000);
  axisGr->Draw("Ap");
  nrGamma->Draw("P");
  nrZmm->Draw("P");
  nrZee->Draw("P");
  TLegend* leg = new TLegend(0.52,0.77,0.88,0.9);
  leg->AddEntry( nrZee,"Z (#rightarrow e^{+}e^{-}) + jet","P");
  leg->AddEntry( nrZmm,"Z (#rightarrow#mu^{+}#mu^{-}) + jet","P");
  leg->AddEntry( nrGamma,"#gamma + jet","P");
  leg->SetMargin(0.2);
  leg->SetFillColor(0);
  leg->Draw();
  TLatex* CMS = new  TLatex();
  CMS->SetTextAlign(12);
  CMS->SetTextSize(0.04);
  CMS->SetNDC();
  CMS->DrawLatex(0.5, 0.32, "CMS Preliminary");
  //  CMS->DrawLatex(0.5, 0.22, "#intLdt = 100 pb^{-1}");
  CMS->DrawLatex(0.52, 0.27, "MC Statistics");
  gPad->SetLogx();
  gPad->SetGridy();
  NormRes->SaveAs("Fig6b.eps");
  NormRes->SaveAs("Fig6b.gif");
  NormRes->SaveAs("Fig6b.root");
}





///////////////////////////////////////////////////////////////////////
double FitUncertainty(TF1* f, TMatrixD* COV, double x)
{
  int dim = COV->GetNrows(),N = f->GetNumberFreeParameters(), npar=5;
  double df,sum,y;
  double PartialDerivative[10],Parameter[10];

  if (dim != npar || N != npar)
    {
      cout<<"ERROR: wrong number of parameters !!!!"<<endl;
      return(-1);
    }  
  for(int i=0;i<npar;i++) Parameter[i] = f->GetParameter(i);
  double z = pow(log10(x),Parameter[2]);  
  PartialDerivative[0] = 1.;
  PartialDerivative[1] = -1./(z+Parameter[3]);
  PartialDerivative[2] = PartialDerivative[3]*log(log10(x))*z;  
  PartialDerivative[3] = Parameter[1]/pow(z+Parameter[3],2);
  PartialDerivative[4] = 1./x;

  sum = 0.;
  for(int i=0;i<npar;i++)
    for(int j=0;j<npar;j++)
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

TH1F* Uncertainty(const TF1& f, TMatrixD& COV, char* name) {
  
  TH1F* h = new TH1F( name, "", 80000, 20, 820);
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
  h->GetXaxis()->SetTitle("p_{T}^{Z/#gamma} (GeV/c)");
  h->GetYaxis()->SetTitle("Response");  
  h->SetLineColor(kGray+2);
  h->SetFillColor(kGray+2);
  h->SetMarkerColor(kGray+2);  

  return h;  
}
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
TGraphErrors* ConvertFromAsymGraph(const TGraphAsymmErrors& g) {

  vector<double> x;
  vector<double> y;
  vector<double> xerr;
  vector<double> yerr;

  for(int i=0; i<g.GetN(); ++i) {
    x.push_back( g.GetX()[i] );
    xerr.push_back( g.GetEXhigh()[i] );
    y.push_back( g.GetY()[i] );
    yerr.push_back( g.GetEYhigh()[i] );
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
TGraphErrors* GetNormResidual(const TGraphErrors& g, const TF1& f) {

  vector<double> x;
  vector<double> y;
  vector<double> xerr;
  vector<double> yerr;

  double res, error;

  for(int i=0; i<g.GetN(); ++i) {
    double xi = g.GetX()[i];
    x.push_back( xi );
    xerr.push_back( g.GetEX()[i] );

    res = 0.0;
    error = 0.0;

    if(f.Eval(xi) != 0.0) {
      res = g.GetY()[i]/f.Eval(xi) - 1.0;
      error =  g.GetEY()[i] / f.Eval(xi);
    }
    y.push_back( res );
    yerr.push_back( error );
  }

  TGraphErrors* new_gr = new TGraphErrors(x.size(),
				     &(x[0]),&(y[0]),
                                     &(xerr[0]),&(yerr[0]));
  new_gr->SetName(g.GetName());
  new_gr->SetTitle("");
  return new_gr;  
}
