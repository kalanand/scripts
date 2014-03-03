// Makes a plot with a likelihood cut:

RooPlot * plot(RooRealVar * var,
	       RooAbsData * dat, 
	       RooAbsPdf * pdf, 
	       const char * logLCut = "nll<-4") {

  // always use the signal likelihood:
  RooAbsPdf * sigPdf = pdfOnResDK.sigGoodD0()->getPdf();
  
  // create a projection of the PDF over the plotting variable:
  RooAbsReal * pdfProj = 
    sigPdf->createProjection(RooArgSet(*Deltae, *d0mass, *nnout, *bknnout), 
			     *var);

  // create the likelihood and add it to the data set as a column:  
  RooFormulaVar nllFunc("nll", "-log(likelihood)", "-log(@0)", *pdfProj);
  RooRealVar * nll = ((RooDataSet*)data)->addColumn(nllFunc);

  // plot the likelihood:
  logLframe = nll->frame(-20, 15);
  data->plotOn(logLframe);

  // apply the log likelihood cut on the data:
  RooDataSet * sliceData = (RooDataSet*)(((RooDataSet*)data)->reduce(logLCut));
  RooPlot * result = var->frame();
  sliceData->plotOn(result);
  pdf->plotOn(result, ProjWData(*sliceData));
}
