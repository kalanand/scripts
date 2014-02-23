// DKbgdStudy() does a full fit toy study, comparing the results of
// the fit before and after adding maximally-biasing DK_bgd sample
// events the fit sample. Naximally biasing means that these are
// events with the good D Dalitz PDF but the wrong-charge kaon. 

RooDataSet * dataNoBias = 0;
TH2 * hist0 = 0;
TH2 * hist1 = 0;
TH2 * hist2 = 0;

void DKbgdStudy(const char * outFileName = "DKbgdStudy.dat", 
		Bool_t noOtherBgd = kFALSE, 
		Bool_t sigOnly = kFALSE,
		int nBiasEvts = 4,
		const char * extraParFileName = 0) {

  // Make sure nBiasEvts is even:
  if ((nBiasEvts / 2) * 2 != nBiasEvts) {
    cerr <"3rd argument to DKbgdStudy is not even. Reducing by 1." << endl;
    nBiasEvts -= 1;
  }

  // This is the # of events to be generated of each charge:
  const int halfNBiasEvts = nBiasEvts / 2;

  ofstream outFile(outFileName);

  // Reread params to recover from previous run and set the yields to
  // the right values:
  pdfOnResDK.parameters().readFromFile("../BToDKTo3piK/params/all.par");
  if (0 != extraParFileName) {
    // read additional parameters:
    cout << "Reading par file " << extraParFileName << endl;
    pdfOnResDK.parameters().readFromFile(extraParFileName);

    // Then need to recalculate interference normalization terms for the B PDF
    // and reevaluate the polar coordinates:
    dalitzHolderN.sigGoodD0Type().recalcX0();
    dalitzHolderN.sigGoodD0Type().setPolarCoords(0.077, 0.064);
    dalitzHolderP.sigGoodD0Type().recalcX0();
    dalitzHolderP.sigGoodD0Type().setPolarCoords(-0.129, 0.019);
  }

  // Set the nsig and asym from the CP parameters:
  pdfOnResDK.setNsigAsymFromXY();
  
  if (noOtherBgd) {
    // Remove non-DKbad bgd:
    ((RooRealVar*)(pdfOnResDK.parameters().
		  find("pdfOnResDK.DpiGoodD0NumEvts")))->setVal(0);

    ((RooRealVar*)(pdfOnResDK.parameters().
		  find("pdfOnResDK.qqBadD0NumEvts")))->setVal(0);

    ((RooRealVar*)(pdfOnResDK.parameters().
		  find("pdfOnResDK.totBBNumEvts")))->setVal(0);
  }

  // Get the polar coords:
  RooRealVar* thetaMinus = 
    (RooRealVar*)(pdfOnResDK.parameters().
		  find("dalitzHolderN.sigGoodD0.theta"));
  
  RooRealVar* rhoMinus = 
    (RooRealVar*)(pdfOnResDK.parameters().
		  find("dalitzHolderN.sigGoodD0.rho"));
  
  RooRealVar* thetaPlus = 
    (RooRealVar*)(pdfOnResDK.parameters().
		  find("dalitzHolderP.sigGoodD0.theta"));
  
  RooRealVar* rhoPlus = 
    (RooRealVar*)(pdfOnResDK.parameters().
		  find("dalitzHolderP.sigGoodD0.rho")); 

  // and their initial, true values:
  const double rhoM = rhoMinus->getVal();
  const double thetaM = thetaMinus->getVal();
  const double rhoP = rhoPlus->getVal();
  const double thetaP = thetaPlus->getVal();

  // Before generating the sample, reduce the DKbgd yield
  // by nBiasEvts, so that these events can be handled separately:

  // First, compute the new DKbgdFrac so as to remove nBiasEvts DK_bgd events: 
  RooRealVar * sigBadD0Frac = (RooRealVar*)
    (pdfOnResDK.parameters().find("pdfOnResDK.sigBadD0Frac"));

  const double oldSigBadD0Frac = sigBadD0Frac->getVal();
  const int nSigGoodD0 = pdfOnResDK.numEvt(BdkEvtTypes::SIG_GOOD_D)->getVal();
  const int nSigBadD0 = oldSigBadD0Frac * nSigGoodD0;
    
  double newSigBadD0Frac = (nSigBadD0 - nBiasEvts) / nSigGoodD0;

  cout << "Initial params:"
       << " oldSigBadD0Frac=" << oldSigBadD0Frac
       << " nSigGoodD0=" << nSigGoodD0
       << " nSigBadD0=" << nSigBadD0 << endl;

  sigBadD0Frac->setVal((double)(nSigBadD0 - nBiasEvts) / nSigGoodD0);

  cout << "To reduce nSigBadD0 by " << nBiasEvts
       << " events, setting sigBadD0Frac to " << sigBadD0Frac->getVal()
       << endl;
    
  BdkPdfAbsBase * pdfBadP = &(pdfOnResDK.sigBadD0P());
  BdkPdfAbsBase * pdfBadN = &(pdfOnResDK.sigBadD0N());

  if (sigOnly) {
    cout << "Generating signal only." << endl;
    
    ((RooRealVar*)(pdfOnResDK.parameters().
		  find("pdfOnResDK.sigBadD0Frac")))->setVal(0);

    pdfBadP = &(pdfOnResDK.sigGoodD0P());
    pdfBadN = &(pdfOnResDK.sigGoodD0N());
  }

  // Given this reduced # of DKbad events, generate the sample:
  cout << "Generating sample without bias DKbgd events" << endl;
  dataNoBias = pdfOnResDK.generate();
  
  // Now generate the missing DKbgd events with their usual PDF (for
  // the first fit):
  cout << "Generating " << nBiasEvts << " DKbgd events." << endl;
  RooDataSet * origBadDataP = pdfBadP->generate(halfNBiasEvts);
  RooDataSet * origBadDataN = pdfBadN->generate(halfNBiasEvts);
  
  // Create a new data set that contains the usual + bias events:
  data = new RooDataSet(*dataNoBias);
  RooAbsData * dataNoBiasP = dataNoBias->reduce("Hdtrkchge>0");
  hist0 = ((RooDataSet*)dataNoBiasP)->createHistogram(*m12, *m13);

  if (0 != origBadDataP) data->append(*origBadDataP);
  if (0 != origBadDataN) data->append(*origBadDataN);

  RooAbsData * dataP = data->reduce("Hdtrkchge>0");
  hist1 = ((RooDataSet*)dataP)->createHistogram(*m12, *m13);

  // Now perform the first fit on a normal, unbiased sample. We are
  // doing the Dalitz plot fit only, not doing a yields fit, so the
  // yields parameters are set to their true values.
  cout << "Doing first fit with unbiased sample" << endl;
  pdfOnResDK.fit(*data, 0.0);
  
  outFile << pdfOnResDK.xyFitResult()->status() << " " 
	  << rhoM << " " 
	  << thetaM << " " 
	  << rhoP << " " 
	  << thetaP << " " 
	  << rhoMinus->getVal() << " " << rhoMinus->getError() << " " 
	  << thetaMinus->getVal() << " " << thetaMinus->getError() << " "
	  << rhoPlus->getVal() << " " << rhoPlus->getError() << " " 
	  << thetaPlus->getVal() << " " << thetaPlus->getError() << " ";
    
  cout << pdfOnResDK.xyFitResult()->status() << " " 
       << rhoM << " " 
       << thetaM << " " 
       << rhoP << " " 
       << thetaP << " " 
       << rhoMinus->getVal() << " " << rhoMinus->getError() << " " 
       << thetaMinus->getVal() << " " << thetaMinus->getError() << " "
       << rhoPlus->getVal() << " " << rhoPlus->getError() << " " 
       << thetaPlus->getVal() << " " << thetaPlus->getError() << " ";

  /*    
  // Now vary the Dalitz variables of the added DKbgd events to be like
  // those of signal but with the wrong kaon charge, so as to simulate 
  // a good-D-bad-K reconstruction:
  cout << "Changing these events' Dalitz vars to wrong-sign signal" << endl;
  RooDataSet * newGoodDataP = pdfOnResDK.sigGoodD0P().generate(halfNBiasEvts);
  RooDataSet * newGoodDataN = pdfOnResDK.sigGoodD0N().generate(halfNBiasEvts);
  
  for (int e = 0; e < halfNBiasEvts; ++e) {
    // the original events and their Dalitz variables:
    RooArgSet * origBadEvtN = origBadDataN->get(e);
    RooArgSet * origBadEvtP = origBadDataP->get(e);
    
    RooRealVar * m12BadN = (RooRealVar*)(origBadEvtN->find(m12->GetName()));
    RooRealVar * m12BadP = (RooRealVar*)(origBadEvtP->find(m12->GetName()));
    
    // The replacement values:    
    RooArgSet * newGoodEvtN = newGoodDataN->get(e);
    RooArgSet * newGoodEvtP = newGoodDataP->get(e);
    
    RooRealVar * m12GoodN = (RooRealVar*)(newGoodEvtN->find(m12->GetName()));
    RooRealVar * m12GoodP = (RooRealVar*)(newGoodEvtP->find(m12->GetName()));
    
    // Do the replacement, flipping charges:
    m12BadN->setVal(m12GoodP->getVal());
    m12BadP->setVal(m12GoodN->getVal());
  }

  // We changed the events in origBadData*. Now copy them to data:
  data = new RooDataSet(*dataNoBias);
  if (0 != origBadDataP) data->append(*origBadDataP);
  if (0 != origBadDataN) data->append(*origBadDataN);
  */

  // Since there is a problem with RooDataSet::get(), we create a new
  // RooProdPdf for these events and regenerate them from scratch.
  // Note that we flip the sign for the Dalitz PDF:

  // The variables to generate:
  RooArgSet genSet(*Deltae, *qprime, *m12, *m13);

  // Setup the N pdf and generate its events:
  RooArgList listN(*(pdfOnResDK.sigBadD0N().getDeltaEPdf()->getPdf()),
		   *(pdfOnResDK.sigBadD0N().getNnContPdf()->getPdf()),
		   *(pdfOnResDK.sigBadD0P().getDalitzPdf()->getPdf()));

  RooProdPdf pdfN("pdfN", "pdfN", listN);

  RooDataSet * newDataN = pdfN->generate(genSet, halfNBiasEvts);
  Hdtrkchge->setIndex(-1);
  newDataN->addColumn(*Hdtrkchge);

  // Setup the P pdf and generate its events:
  RooArgList listP(*(pdfOnResDK.sigBadD0P().getDeltaEPdf()->getPdf()),
		   *(pdfOnResDK.sigBadD0P().getNnContPdf()->getPdf()),
		   *(pdfOnResDK.sigBadD0N().getDalitzPdf()->getPdf()));

  RooProdPdf pdfP("pdfP", "pdfP", listP);

  RooDataSet * newDataP = pdfP->generate(genSet, halfNBiasEvts);
  Hdtrkchge->setIndex(1);
  newDataP->addColumn(*Hdtrkchge);

  // put it all into the data:
  data = new RooDataSet(*dataNoBias);

  if (0 != newDataP) data->append(*newDataP);
  if (0 != newDataN) data->append(*newDataN);

  // What does it look like for one charge:
  dataP = data->reduce("Hdtrkchge>0");
  hist2 = ((RooDataSet*)dataP)->createHistogram(*m12, *m13);

  // Refit:
  cout << "Doing second fit with bias" << endl;
  pdfOnResDK.fit(*data, 0.0);
  
  outFile << pdfOnResDK.xyFitResult()->status() << " " 
	  << rhoMinus->getVal() << " " << rhoMinus->getError() << " " 
	  << thetaMinus->getVal() << " " << thetaMinus->getError() << " "
	  << rhoPlus->getVal() << " " << rhoPlus->getError() << " " 
	  << thetaPlus->getVal() << " " << thetaPlus->getError() << endl;

  cout << pdfOnResDK.xyFitResult()->status() << " " 
       << rhoMinus->getVal() << " " << rhoMinus->getError() << " " 
       << thetaMinus->getVal() << " " << thetaMinus->getError() << " "
       << rhoPlus->getVal() << " " << rhoPlus->getError() << " " 
       << thetaPlus->getVal() << " " << thetaPlus->getError() << endl;

  TCanvas * can = new TCanvas("can", "can", 1200, 400);
  can->Divide(3,1);
  can->cd(1);
  hist0->Draw();
  can->cd(2);
  hist1->Draw();
  can->cd(3);
  hist2->Draw();
}

