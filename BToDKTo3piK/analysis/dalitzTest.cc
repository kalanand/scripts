// $Id: dalitzTest.cc,v 1.10 2006/04/28 01:56:00 fwinkl Exp $
// Script to test the daliz integration and boundaries


void dalitzTest(Bool_t boundaryTest = false, Bool_t invert = kFALSE)
{
  RooNumIntConfig* cfg = RooAbsReal::defaultIntegratorConfig();
  cfg->setEpsAbs(1E-5);
  cfg->setEpsRel(1E-5);
  cfg->method1D().setLabel("RooSegmentedIntegrator1D");

  // RooRealVar m12("m12","",1.5,0,3);
  // RooRealVar m13("m13","",1.5,0,3);

  /*
  //RooRealVar a0("a0","",6.25247e+00);
  RooRealVar a0("a0","",0);
  //  RooRealVar a1("a1","",-1.33799e+01);
  RooRealVar a1("a1","",0);
  RooRealVar a2("a2","",-1.54353e+01);
  RooRealVar a3("a3","",4.42504e+00);
  RooRealVar a4("a4","",7.61593e+00);
  RooRealVar a5("a5","",8.72822e+00);
  RooRealVar a6("a6","",-4.41505e-01);
  RooRealVar a7("a7","",-1.23653e+00);
  RooRealVar a8("a8","",-1.10855e+00);
  //  RooRealVar a9("a9","",-2.52116e+00);
  RooRealVar a9("a9","",0);
  */
  
  RooRealVar a0("a0","",1);
  RooRealVar a1("a1","",0);
  RooRealVar a2("a2","",0);
  RooRealVar a3("a3","",0);
  RooRealVar a4("a4","",0);
  RooRealVar a5("a5","",0);
  RooRealVar a6("a6","",0);
  RooRealVar a7("a7","",0);
  RooRealVar a8("a8","",0);
  RooRealVar a9("a9","",0);
  
  RooRealVar *a[10] = {&a0,&a1,&a2,&a3,&a4,&a5,&a6,&a7,&a8,&a9};

  
  // Need three pdfs here since RooFit caches the normalization
  BdkPdf2DpolyDalitz pdf("pdf","",*s12,*s13,1,a);
  BdkPdf2DpolyDalitz pdf2("pdf","",*s12,*s13,1,a);
  BdkPdf2DpolyDalitz pdf3("pdf","",*s12,*s13,1,a);
  
  RooAbsPdf::verboseEval(0);

  /*
  for (int i = 0;i<100;i++) {
    *s13.setVal(3.0/100*i);    
    cout << "Value of PDF at "<< *s13.getVal()<<" "<<*s12.getVal()<<" "<< pdf.getPdf()->getVal()<<endl;
  }
  
  RooDataSet *data = pdf.generate(500);
  pdf.getPdf()->forceNumInt();
  pdf.getPdf()->verboseEval(0);
  RooPlot *plot = *s12.frame();
  data->plotOn(plot);
  pdf.getPdf()->plotOn(plot);
  plot->Draw();
  */

  // Pure numerical integration
  pdf2.getPdf()->forceNumInt();
  pdf2.getPdf()->setIntegratorConfig(*cfg);
  Double_t num1 = pdf2.getPdf()->getNorm(&RooArgSet(*s12));
  Double_t num2 = pdf2.getPdf()->getNorm(&RooArgSet(*s13));
  Double_t num3 = pdf2.getPdf()->getNorm(&RooArgSet(*s12,*s13));
  
  cout << "numerical int over m12  "<<num1<<endl;
  cout << "numerical int over m13  "<<num2<<endl;
  cout << "numerical int over both "<<num3<<endl;

  // Pure analytical integration
  ((Bdk2DpolyDalitz*)pdf3.getPdf())->customInt(1);
  Double_t ana1 = pdf3.getPdf()->getNorm(&RooArgSet(*s12,*s13));
  cout << "analytical int over both " << ana1 <<endl;

  // Hybrid integrations
  RooAbsPdf::verboseEval(0);
  //pdf.setM23VetoMass(0.4,0.401);
  Double_t custom1 = pdf.getPdf()->getNorm(&RooArgSet(*s12));
  Double_t custom2 = pdf.getPdf()->getNorm(&RooArgSet(*s13));
  Double_t custom3 = pdf.getPdf()->getNorm(&RooArgSet(*s12,*s13));

  cout << "custom int over m12 " << custom1 <<endl;
  cout << "custom int over m13 " << custom2 <<endl;
  cout << "custom int over both " << custom3 <<endl;

  //((BdkDalitzBase*)pdf.getPdf())->setMasses(2,0.1,0.2,0.4);

  if (boundaryTest) {
    const int POINTS = 3000;

    RooAbsPdf::verboseEval(0);
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c","dalitzTest",600,600);
    TGraph* g2 = new TGraph(POINTS);
   
    // This tests the drawBoundary method
    TGraph *g = ((BdkDalitzBase*)pdf.getPdf())->drawBoundary(1000); 
    if (invert) {
      g->SetLineColor(kBlue);
      g->SetLineWidth(2);
    }
    else {
      g->SetLineColor(kRed);
      g->SetLineWidth(1);
    }

    int p = 0;
    // This tests the inDalitz method
    for (Double_t x=0;x<=3;x+=(3.0/POINTS)) {
      for (Double_t y=0;y<=3;y+=(3.0/POINTS)) {
        if (invert) {
          if (((BdkDalitzBase*)pdf.getPdf())->inM23Veto(x,y) &&
              ((BdkDalitzBase*)pdf.getPdf())->inDalitzBounds(x,y))
            g2->SetPoint(p++,x,y);
        }
        else {
          if (((BdkDalitzBase*)pdf.getPdf())->inDalitz(x,y))
            g2->SetPoint(p++,x,y);
        }
      }
    }
    
    g2->SetMarkerColor(kBlue);
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(0.1);

    g->Draw("ac");
    g2->Draw("p same");

    g->GetXaxis()->SetTitle(m12->GetTitle());
    g->GetYaxis()->SetTitle(m13->GetTitle());
    g->SetTitle("");
    c->Update();
    c->SaveAs("dalitzTest.eps");
  }
}
