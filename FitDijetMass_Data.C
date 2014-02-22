
/*******************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: Presently in the user code
 *
 *
 * Authors:
 *
 *   Kalanand Mishra, Fermilab - kalanand@fnal.gov
 *
 * Description:
 *   Plots dijet mass distributions in data and fit to QCD lineshape 
 *
 * Implementation details:
 *   Uses RooFit and RooStat classes.
 *  
 * History:
 *   
 *
 * Copyright (C) 2010 FNAL 
 ********************************************************************/


#include <memory>
#include "QCDFitUncertainty.h"

//  _     _ _                          
// | |   (_) |__  _ __ __ _ _ __ _   _ 
// | |   | | '_ \| '__/ _` | '__| | | |
// | |___| | |_) | | | (_| | |  | |_| |
// |_____|_|_.__/|_|  \__,_|_|   \__, |
//                               |___/ 

gSystem->Load("libRooFit") ;
using namespace RooFit;



//  __  __       _         _____                 _   _             
// |  \/  | __ _(_)_ __   |  ___|   _ _ __   ___| |_(_) ___  _ __  
// | |\/| |/ _` | | '_ \  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \ 
// | |  | | (_| | | | | | |  _|| |_| | | | | (__| |_| | (_) | | | |
// |_|  |_|\__,_|_|_| |_| |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|
// 

void FitDijetMass_Data() {

  
  TFile *inf  = new TFile("MassResults_ak7calo.root");
  TH1F *hCorMassDen     = (TH1F*) inf->Get("DiJetMass");
  hCorMassDen->SetXTitle("Corrected Dijet Mass (GeV)");
  hCorMassDen->SetYTitle("Events/GeV");
  hCorMassDen->GetYaxis()->SetTitleOffset(1.5);
  hCorMassDen->SetMarkerStyle(20);
  hCorMassDen->GetXaxis()->SetRangeUser(120.,900.);



  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPadRightMargin(0.08);
  tdrStyle->SetLegendBorderSize(0);
  gStyle->SetOptFit(1111);
  tdrStyle->SetOptStat(0); 

  /*
  TCanvas* c2 = new TCanvas("c2","DijetMass", 500, 500);
  /////// perform 4 parameters fit
  TF1 *func = new TF1("func", "[0]*((1-x/7000.+[3]*(x/7000)^2)^[1])/(x^[2])", 
  100., 1000.);
  func->SetParameter(0, 1.0e+08);
  func->SetParameter(1, -1.23);
  func->SetParameter(2, 4.13);
  func->SetParameter(3, 1.0);

  func->SetLineColor(4);
  func->SetLineWidth(3);

  TVirtualFitter::SetMaxIterations( 10000 );
  TVirtualFitter *fitter;
  TMatrixDSym* cov_matrix;

  int fitStatus = hCorMassDen->Fit("func","LLI","",130.0, 800.0); // QCD fit
 
  TH1F *hFitUncertainty = hCorMassDen->Clone("hFitUncertainty");
  hFitUncertainty->SetLineColor(5);
  hFitUncertainty->SetFillColor(5);
  hFitUncertainty->SetMarkerColor(5);

  if (fitStatus == 0) {
    fitter = TVirtualFitter::GetFitter();
    double* m_elements = fitter->GetCovarianceMatrix();
    cov_matrix = new TMatrixDSym( func->GetNumberFreeParameters(),m_elements);
    cov_matrix->Print();
    double x, y, e;

    for(int i=0;i<hFitUncertainty->GetNbinsX();i++)
      {
	x = hFitUncertainty->GetBinCenter(i+1);
	y = func->Eval(x);
	e = QCDFitUncertainty( func, *cov_matrix, x);
	hFitUncertainty->SetBinContent(i+1,y);
	hFitUncertainty->SetBinError(i+1,e);
      }
  }

  hCorMassDen->Draw("ep");
  gPad->Update();
  TPaveStats *st = (TPaveStats*)hCorMassDen->FindObject("stats");
  st->SetName("stats1");
  st->SetX1NDC(0.3); //new x start position
  st->SetX2NDC(0.6); //new x end position
  st->SetTextColor(4);
  hCorMassDen->GetListOfFunctions()->Add(st);


  /////// perform 2 parameters fit
  TF1 *func2 = new TF1("func2", "[0]*(1-x/7000.)/(x^[1])", 100., 1000.);
  func2->SetParameter(0, 10000.);
  func2->SetParameter(1, 5.0);
  func2->SetLineWidth(3);

  fitStatus = hCorMassDen->Fit("func2","LLI","",130.0, 800.0); // QCD fit

  TH1F *hFitUncertainty2 = hCorMassDen->Clone("hFitUncertainty2");
  hFitUncertainty2->SetLineColor(kGray);
  hFitUncertainty2->SetFillColor(kGray);
  hFitUncertainty2->SetMarkerColor(kGray);

  if (fitStatus == 0) {
    fitter = TVirtualFitter::GetFitter();
    double* m_elements = fitter->GetCovarianceMatrix();
    cov_matrix = new TMatrixDSym( func2->GetNumberFreeParameters(),m_elements);
    cov_matrix->Print();
    double x, y, e;

    for(int i=0;i<hFitUncertainty2->GetNbinsX();i++)
      {
	x = hFitUncertainty2->GetBinCenter(i+1);
	y = func2->Eval(x);
	e = QCDFitUncertainty( func2, *cov_matrix, x);
	hFitUncertainty2->SetBinContent(i+1,y);
	hFitUncertainty2->SetBinError(i+1,e);
      }
  }

  hFitUncertainty->Draw("E3 same");
  hCorMassDen->Draw("ep sames");
  hFitUncertainty2->Draw("E3 same");
  hCorMassDen->Draw("ep sames");
  func2->Draw("same");
  c2->SetLogy(1);

  */


  

  TH1F *hCorMass     = hCorMassDen->Clone("hCorMass");


  for(int i=0; i<hCorMass->GetNbinsX(); i++){
    hCorMass->SetBinContent(i+1, hCorMassDen->GetBinContent(i+1)*hCorMassDen->GetBinWidth(i+1));
    hCorMass->SetBinError(i+1, hCorMassDen->GetBinError(i+1) * hCorMassDen->GetBinWidth(i+1));
  } 




  // Our observable is the invariant mass
  RooRealVar invMass("invMass", "Corrected dijet mass", 1., 7000.0, "GeV");
  RooDataHist data( "data", "", invMass, hCorMass);

   //////////////////////////////////////////////




   // make QCD model
  RooRealVar p0("p0", "# events", 100000., -10., 1.0e+8);
  RooRealVar p1("p1","p1", 0.) ;  
  RooRealVar p2("p2","p2", 4.8, -100., 100.) ; 
  RooRealVar p3("p3","p3", 0.) ; 



//    // define QCD line shape
  RooGenericPdf qcdModel("qcdModel", "pow(1-@0/7000.+@3*(@0/7000.)*(@0/7000.),@1)/pow(@0,@2)",RooArgList(invMass,p1,p2,p3)); 



   // full model
  //  RooAddPdf model("model","qcd",RooArgList(qcdModel), RooArgList(p0)); 
  invMass.setRange("fitRange",137.,1000.);
  RooExtendPdf model("model",  "", qcdModel,   p0, "fitRange");



   //plot sig candidates, full model, and individual componenets 

   //   __ _ _    
   //  / _(_) |_  
   // | |_| | __| 
   // |  _| | |_  
   // |_| |_|\__| 


 // Important: fit integrating f(x) over ranges defined by X errors, rather
  // than taking point at center of bin

   RooFitResult* fit = model.fitTo(data, Minos(kFALSE), Hesse(kFALSE),Extended(kTRUE),
   SumW2Error(kFALSE),Save(kTRUE), Range("fitRange"), Integrate(kTRUE) );



   // to perform chi^2 minimization fit instead
   //    RooFitResult* fit = model.chi2FitTo(data, Extended(kTRUE), 
   // 				       Save(),Range(50.,526.),Integrate(kTRUE) );

   fit->Print();
 

   const RooArgSet* normSet = new RooArgSet(invMass);
//    cout << "pdf norm = " << model.expectedEvents(normSet) << endl;
//    cout << "data norm = " << hCorMass->Integral() << endl;

   int firstBin = hCorMass->FindFirstBinAbove(0.0);
   invMass.setVal(hCorMass->GetBinLowEdge(firstBin));
   double a = model.getVal(normSet);
   invMass.setVal(hCorMass->GetBinLowEdge(firstBin+1));
   double b = model.getVal(normSet); 
   // invMass.setRange("firstBinRange",a,b);

   // double relNorm = 0.27;
   double relNorm = 0.5*(a-b) * model.expectedEvents(normSet) * 
      hCorMass->GetBinWidth(firstBin) / hCorMass->Integral(firstBin, firstBin);

//    cout << "pdf first bin = " << pdfValFirstBin << endl;
//    cout << "data first bin = " << hCorMass->Integral(firstBin, firstBin) << endl;
   // invMass.setRange("fitRange");

   //plot data 
   TCanvas* cdataNull = new TCanvas("cdataNull","fit to dijet mass",500,500);
   RooPlot* frame1 = invMass.frame(0., 1000.) ; 
   data.plotOn(frame1, DataError(RooAbsData::SumW2) ) ; 


   model.plotOn(frame1, Normalization(relNorm, RooAbsReal::Relative) );
   model.plotOn(frame1, VisualizeError(*fit, 1),FillColor(kYellow), 
                Normalization(relNorm, RooAbsReal::Relative) ) ; 
   model.plotOn(frame1, Normalization(relNorm, RooAbsReal::Relative) );

   model.paramOn(frame1, &data );
 
   TPaveText* dataPave = (TPaveText*) frame1->findObject("model_paramBox");
   dataPave->SetY1(0.77);
   char evRange[200];
   TPaveText *tbox = new TPaveText(0.55,0.65,0.9,0.75,"BRNDC");
   sprintf(evRange,"%s%lf %s %d","#chi^{2}/#nu = ", frame1->chiSquare(), 
   "for #nu=", fit->floatParsFinal().getSize() );
   tbox->AddText(evRange); 
   tbox->SetBorderSize(1);  
   frame1->addObject(tbox);
   gPad->SetLogy();
   frame1->GetXaxis()->SetRangeUser(80, 1000);
   // frame1->GetYaxis()->SetNoExponent();
   frame1->GetXaxis()->SetNdivisions(505);
   frame1->GetYaxis()->SetRangeUser(1E-1,1E+5);
   frame1->GetYaxis()->SetTitle("Events / bin");
   frame1->GetYaxis()->SetTitleOffset(1.35);
   frame1->SetTitle("fit to data with QCD lineshape");
   frame1->Draw() ;



    
    // S h o w   r e s i d u a l   a n d   p u l l   d i s t s
    // -------------------------------------------------------
    
   //// Construct a histogram with the residuals of the data w.r.t. the curve
   RooHist* hresid = frame1->residHist() ;
   // Create a new frame to draw the residual distribution and add the distribution to the frame
   RooPlot* frame2 = invMass.frame() ;
   frame2->GetYaxis()->SetTitle("Data - Fit");
   frame2->GetXaxis()->SetRangeUser(80, 1000);
   frame2->addPlotable(hresid,"P") ;


    
   ///// Construct a histogram with the pulls of the data w.r.t the curve
   RooHist* hpull = frame1->pullHist() ;   
   //// Create a new frame to draw the pull distribution and add the distribution to the frame
   RooPlot* frame3 = invMass.frame(Title("Pull Distribution")) ;
   frame3->GetYaxis()->SetTitle("(Data - Fit) / error");
   frame3->GetXaxis()->SetRangeUser(80, 1000);
   frame3->addPlotable(hpull,"P") ;


   TCanvas* cResidual = new TCanvas("cResidual","Residual Distribution",1000,500);
   cResidual->Divide(2) ;
   cResidual->cd(1) ; 
   gPad->SetLeftMargin(0.15) ; 
   // frame1->GetYaxis()->SetTitleOffset(1.6) ; 
   frame2->Draw() ;
   cResidual->cd(2) ; 
   gPad->SetLeftMargin(0.15) ; 
   // frame1->GetYaxis()->SetTitleOffset(1.6) ; 
   frame3->Draw() ;
  



  

}






