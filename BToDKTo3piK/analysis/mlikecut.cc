void mlikecut() {
                                                                                                                     
  RooAbsPdf * pdf = pdfOnResDK.sigGoodD0()->getPdf();
  RooAbsPdf * totPdf = pdfOnResDK.getPdf();                                                                                                                   
  // create a projection of the PDF over mes (the plotting variable):
  RooAbsReal * pdfProj = pdf->createProjection(RooArgSet(*Deltae, *d0mass, *nnout, *bknnout), *mes);
                                                                                                                     
  // create the likelihood and add it to the data set as a column:
  RooFormulaVar nllFunc("nllFunc", "-log(likelihood)", "-log(@0)", *pdfProj);
  RooRealVar * nll = ((RooDataSet*)data)->addColumn(nllFunc);
                                                                                                                     
  // plot the likelihood:
  //RooPlot * pframe = nll->frame(-20, 15);
  //data->plotOn(pframe);
  //pframe->Draw();                                                                                                                     
  // apply the log likelihood cut on the data:
  RooDataSet * sliceData =  (RooDataSet*)(((RooDataSet*)data)->reduce("nllFunc<-5.5"));

  RooArgSet chargeSet(*Hdtrkchge);
  RooDataSet redData("redData", "reduced Data for plotting", chargeSet);  
  Hdtrkchge->setIndex(-1);
  redData.add(chargeSet);
  Hdtrkchge->setIndex(1);
  redData.add(chargeSet);

  sliceData->plotOn(mesFrame);
  totPdf->plotOn(mesFrame,ProjWData(*sliceData));
  mesFrame->getAttLine()->SetLineColor(kRed);
  RooArgSet args(*(pdfOnResDK.qqBadD0()->getPdf()));
  args.add(*(pdfOnResDK.qqGoodD0()->getPdf()));
  totPdf->plotOn(mesFrame,Components(args),ProjWData(*sliceData));
  mesFrame->getAttLine()->SetLineColor(kBlack); 
  args.add(*(pdfOnResDK.BBBadD0()->getPdf()));
  args.add(*(pdfOnResDK.BBGoodD0()->getPdf()));
  totPdf->plotOn(mesFrame,Components(args),ProjWData(*sliceData));
  mesFrame->getAttLine()->SetLineColor(kGreen); 
  args.add(*(pdfOnResDK.DpiGoodD0()->getPdf()));
  args.add(*(pdfOnResDK.DpiBadD0()->getPdf()));
  args.add(*(pdfOnResDK.DPiPi()->getPdf()));
  args.add(*(pdfOnResDK.charmless()->getPdf()));
  totPdf->plotOn(mesFrame,Components(args),ProjWData(*sliceData));
  mesFrame->getAttLine()->SetLineColor(kBlue);
 TCanvas * cant = new TCanvas("cant", "cant", 600, 400);
  mesFrame->Draw();
  cant->Print("mlikecut.eps");
}

