// Make the standard plots of the fit vars.
// The canvas is saved as plotFileName (see plot.hh):

#include "../BToDKTo3piK/globals/vars.hh"

void stdPlot(RooAbsData * dat = data, Bool_t plotComps = kTRUE) {

  if (kFALSE == doPlot || 0 == dat) {
    return;
  }
  
  // Plot: 
  RooArgSet chargeSet(*Hdtrkchge);
  RooDataSet catData("catData", "reduced Data for plotting", chargeSet);
  Hdtrkchge->setIndex(-1);
  catData.add(chargeSet);
  Hdtrkchge->setIndex(1);
  catData.add(chargeSet);
  
  setupPlots();

  plot1dIntCfg.setEpsAbs(1e-4);
  plot1dIntCfg.setEpsRel(1e-4);
  pdfOnResDK.getPdf()->setIntegratorConfig(plot1dIntCfg);

  // plot data:
  dat->plotOn(DeltaeFrame);  
  dat->plotOn(qprimeFrame);
  dat->plotOn(m12Frame);
  dat->plotOn(m13Frame);
  
  // overlay continuum PDF:
  RooArgSet args;
  args.add(*(pdfOnResDK.qqBadD0N().getPdf()));
  args.add(*(pdfOnResDK.qqBadD0P().getPdf()));
  args.add(*(pdfOnResDK.qqGoodD0N().getPdf()));  
  args.add(*(pdfOnResDK.qqGoodD0P().getPdf()));  

  if (plotComps) {
    pdfOnResDK.getPdf()->plotOn(DeltaeFrame, Components(args), ProjWData(catData),
                                LineColor(kBlack));
    pdfOnResDK.getPdf()->plotOn(qprimeFrame, Components(args), ProjWData(catData),
                                LineColor(kBlack));
    pdfOnResDK.getPdf()->plotOn(m12Frame, Components(args), ProjWData(catData),
                                LineColor(kBlack));
    pdfOnResDK.getPdf()->plotOn(m13Frame, Components(args), ProjWData(catData),
                                LineColor(kBlack));
  }

  // overlay total bgd:
  args.add(*(pdfOnResDK.BBBadD0N().getPdf()));  
  args.add(*(pdfOnResDK.BBBadD0P().getPdf()));  
  args.add(*(pdfOnResDK.BBGoodD0N().getPdf()));  
  args.add(*(pdfOnResDK.BBGoodD0P().getPdf()));  
  args.add(*(pdfOnResDK.DKXN().getPdf()));
  args.add(*(pdfOnResDK.DKXP().getPdf()));
  args.add(*(pdfOnResDK.DPiXN().getPdf()));
  args.add(*(pdfOnResDK.DPiXP().getPdf()));
  args.add(*(pdfOnResDK.DpiBadD0N().getPdf()));
  args.add(*(pdfOnResDK.DpiBadD0P().getPdf()));
  args.add(*(pdfOnResDK.DpiGoodD0N().getPdf()));
  args.add(*(pdfOnResDK.DpiGoodD0P().getPdf()));
  args.add(*(pdfOnResDK.sigBadD0N().getPdf()));
  args.add(*(pdfOnResDK.sigBadD0P().getPdf()));

  if (plotComps) {
    pdfOnResDK.getPdf()->plotOn(DeltaeFrame, Components(args), ProjWData(catData),
                                LineColor(kBlue));
    pdfOnResDK.getPdf()->plotOn(qprimeFrame, Components(args), ProjWData(catData),
                                LineColor(kBlue));
    pdfOnResDK.getPdf()->plotOn(m12Frame, Components(args), ProjWData(catData),
                                LineColor(kBlue));
    pdfOnResDK.getPdf()->plotOn(m13Frame, Components(args), ProjWData(catData),
                                LineColor(kBlue));
  }

  // overlay total PDF:
  
  pdfOnResDK.getPdf()->plotOn(DeltaeFrame, ProjWData(catData),
                              LineColor(kRed));
  pdfOnResDK.getPdf()->plotOn(qprimeFrame, ProjWData(catData),
                              LineColor(kRed));
  pdfOnResDK.addPdfPos()->plotOn(m12Frame, ProjWData(catData),
                              LineColor(kRed));

  TCanvas *can = new TCanas("can","stdPlot",1000,1000);
  can->Divide(2,2);

  can->cd(1);
  DeltaeFrame->Draw();

  can->cd(2);
  qprimeFrame->Draw();

  can->cd(3);
  m12Frame->Draw();

  can->cd(4);
  m13Frame->Draw();

  cout << "chi^2 of deltae plot = " << DeltaeFrame->chiSquare() << endl;
  cout << "chi^2 of q' plot = " << qprimeFrame->chiSquare() << endl;
  cout << "chi^2 of m12 plot = " << m12Frame->chiSquare() << endl;
  cout << "chi^2 of m13 plot = " << m13Frame->chiSquare() << endl;

  can->SaveAs("stdPlot.eps");
}



//------------------------------------------------------------------
// plot one variable with a signal likelihood cut:
RooPlot * plotLikeCut(RooRealVar * var, 
		      RooAbsData * dat = data,
		      Bool_t plotNll = kFALSE) {

  if (kFALSE == doPlot || 0 == dat || 0 == var) {
    return;
  }
  RooPlot * varFrame = var->frame(30);
  
  RooAbsPdf * pdf = pdfOnResDK.sigGoodD0()->getPdf();
  RooAbsPdf * totPdf = pdfOnResDK.getPdf();      

  // create a projection of the PDF over var:
  RooArgSet fullSet(*Deltae,*nnout, *bknnout);
  fullSet.remove(*var);

  RooAbsReal * pdfProj = pdf->createProjection(fullSet, *var);
  
  // create the likelihood and add it to the data set as a column:
  RooFormulaVar nllFunc("nllFunc", "-log(likelihood)", "-log(@0)", *pdfProj);
  RooRealVar * nll = ((RooDataSet*)dat)->addColumn(nllFunc);

  if (plotNll) {
    // plot the likelihood:
    RooPlot * pframe = nll->frame(-15, 40);
    dat->plotOn(pframe);
    pframe->Draw();     
    return varFrame;
  }  

  // apply the log likelihood cut on the dat:
  RooDataSet * sliceData = 0; 
  if( var == nnout || var == bknnout  ) {
     sliceData =
         (RooDataSet*)(((RooDataSet*)dat)->reduce("nllFunc<-3.0"));
  } else {
    if(var == d0mass ) { 
        sliceData =
         (RooDataSet*)(((RooDataSet*)dat)->reduce("nllFunc<-3.0"));
    } else {
	sliceData =
         (RooDataSet*)(((RooDataSet*)dat)->reduce("nllFunc<-1.0"));
   } 
  }
  
  
  
  RooArgSet chargeSet(*Hdtrkchge);
  RooDataSet catData("catData", "reduced Data for plotting", chargeSet);  
  Hdtrkchge->setIndex(-1);
  catData.add(chargeSet);
  Hdtrkchge->setIndex(1);
  catData.add(chargeSet);

  // plot data:
  sliceData->plotOn(varFrame);

  // plot total pdf:
  totPdf->plotOn(varFrame,ProjWData(*sliceData));
  varFrame->getAttLine()->SetLineColor(kRed);

  // plot continuum pdf:
  RooArgSet args(*(pdfOnResDK.qqBadD0()->getPdf()));
  args.add(*(pdfOnResDK.qqGoodD0()->getPdf()));
  totPdf->plotOn(varFrame,Components(args),ProjWData(*sliceData));
  varFrame->getAttLine()->SetLineColor(kBlack); 

  // plot BB pdf:
  args.add(*(pdfOnResDK.BBBadD0()->getPdf()));
  args.add(*(pdfOnResDK.BBGoodD0()->getPdf()));
  totPdf->plotOn(varFrame,Components(args),ProjWData(*sliceData));
  varFrame->getAttLine()->SetLineColor(kGreen); 

  // plot peaking BB pdf:
  args.add(*(pdfOnResDK.DpiGoodD0()->getPdf()));
  args.add(*(pdfOnResDK.DpiBadD0()->getPdf()));
  args.add(*(pdfOnResDK.DPiPi()->getPdf()));
  args.add(*(pdfOnResDK.charmless()->getPdf()));
  totPdf->plotOn(varFrame,Components(args),ProjWData(*sliceData));
  varFrame->getAttLine()->SetLineColor(kBlue);

  varFrame->Draw();

  return varFrame;
}

//------------------------------------------------------------------
// plot all vars with a signal likelihood cut:
void stdPlotLikeCut(RooAbsData * dat = data) {
  TCanvas *c1 = new TCanvas("c1","mes",600,600);
  TCanvas *c2 = new TCanvas("c2","Deltae",600,600);
  TCanvas *c3 = new TCanvas("c3","d0mass",600,600);
  TCanvas *c4 = new TCanvas("c4","nnout",600,600);
  TCanvas *c5 = new TCanvas("c5","bknnout",600,600);
  
 //c1->cd();
 // mesFrameLike = plotLikeCut(mes, dat);
  c2->cd();
  DeltaeFrameLike = plotLikeCut(Deltae, dat);
  c3->cd();
  d0massFrameLike = plotLikeCut(d0mass, dat);
  c4->cd();
  qprimeFrameLike = plotLikeCut(nnout, dat);
  c5->cd();
  bknnFrameLike = plotLikeCut(bknnout, dat);
}
  
