{
  using namespace RooFit;


// C r e a t e   e n d   p o i n t   p d f   s h a p e s
  // ------------------------------------------------------

  const double xmin = -3.5;
  const double xmax = 3.5;
  const int nbins = (int) ((xmax - xmin)*2);


  // Observable
  RooRealVar x("x", "x", xmin, xmax);


// Lower end point shape: a Gaussian at 0
  RooRealVar gaus1_mean("gaus1_mean","gaus1_mean",0) ;
  RooRealVar gaus1_sigma("gaus1_sigma","gaus1_sigma", 1.0) ;
  RooGaussian gaus1("gaus1","gaus1",x,gaus1_mean,gaus1_sigma) ;

// Lower end point shape: second Gaussian at 0, different width
  RooRealVar gaus1prime_sigma("gaus1prime_sigma","gaus1prime_sigma", 3.0) ;
  RooGaussian gaus1prime("gaus1prime","gaus1prime",x,gaus1_mean,gaus1prime_sigma) ;


  // Create interpolation variable for lower end point shape
  RooRealVar alpha("alpha","alpha", 0.5, 0., 1.) ;
  x.setBins(10,"cache") ;
  alpha.setBins(10,"cache") ;
  // Finally we have a lower end point shape which is a morphed Gaussian
  RooIntegralMorph g1("g1","lower end point shape",gaus1, gaus1prime,x,alpha) ;
  g1.setCacheAlpha(kTRUE) ;

  // Upper end point shape: a narrow Gaussian at 5
  RooGaussian g2("g2","g2",x, RooFit::RooConst(5),RooFit::RooConst(1.0) ) ;





// C r e a t e   i n t e r p o l a t i n g   p d f 
  // -----------------------------------------------

  // Create interpolation variable
  RooRealVar beta("beta","beta", 0.8, 0., 1.) ;
  beta.setBins(10,"cache") ;

  // Construct interpolating pdf in (x,a) represent g1(x) at a=a_min
  // and g2(x) at a=a_max
  RooIntegralMorph model("model","signal shape after morphing",g1,g2,x,beta) ;
  model.setCacheAlpha(kTRUE) ;


// Indicate to the RooAbsCachedPdf base class that for the filling of the
//  cache the traversal of the x should be in the innermost loop, to minimize
//  recalculation of the one-dimensional internal cache for a fixed value of alpha

//  model.preferredObservableScanOrder( RooArgSet(alpha, beta), RooArgSet(beta, alpha) );




//set the precision for the integral of the model here: 10^-6 should be fine, default is 10^-8.
  RooNumIntConfig* cfg = RooAbsReal::defaultIntegratorConfig();
  cfg->setEpsAbs(1E-5);
  cfg->setEpsRel(1E-5);
  model.setIntegratorConfig(*cfg);



 // F i t   p d f   t o   d a t a s e t   
  // -------------------------------------

  // Generate a toy dataset at alpha = 0.5, beta =0.8
  RooDataSet* toydata = model.generate(x,10000);

  // convert unbinned data to a histogram
  TH1D* hist = new TH1D("hist","", nbins, xmin, xmax);
  toydata->tree()->Draw("x>>hist","","goff");

  // convert the histogram to binned data
  RooDataHist data( "data", "", x, hist);

 // Fit pdf to toy data
  //model.fitTo(data,Verbose(kTRUE)) ;
  RooFitResult* r = model.fitTo(data, Save(true), Minos(false), Hesse(false) );
  r->Print();



  // Plot fitted pdf and data overlaid
  TCanvas* c = new TCanvas("mycanvas","Fit using template morphing",800,600) ;
  RooPlot* frame2 = x.frame( Title(""), Bins(nbins), Range(xmin, xmax)) ;
  data.plotOn(frame2) ;
  model.plotOn(frame2) ; 
  gaus1.plotOn(frame2,LineStyle(kDashed), LineColor(kRed)) ; 
  gaus1prime.plotOn(frame2,LineStyle(kDashed), LineColor(kRed)) ; 
  g1.plotOn(frame2,LineColor(kRed), LineWidth(1)) ; 
  g2.plotOn(frame2,LineColor(kBlack), LineWidth(1)) ; 
  model.paramOn(frame2, &data );
  //  model.setCacheAlpha(kFALSE) ;
  frame2->Draw() ;
  c->SaveAs("TemplateMorphing.gif");
}
