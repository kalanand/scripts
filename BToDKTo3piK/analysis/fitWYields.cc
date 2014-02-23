// Script to test the yields chi^2 term in the NLL

void fitWYields(Bool_t useShape  = kFALSE, 
		Bool_t useYields = kTRUE, 
		Bool_t useBoth   = kFALSE,
		TRandom * rand = 0,  // for fluctuating nsig and asym
		Bool_t toy = kTRUE,
                Int_t nEvents = 0) {

  pdfOnResDK.parameters().readFromFile("../BToDKTo3piK/params/all.par");

  if (toy) {
    pdfOnResDK.setNsigAsymFromXY();  
  }

  if (0 != rand) {
    // fluctuate the signal and asymmetry yields in a simple,
    // uncorrelated way:
    RooRealVar * nsig = pdfOnResDK.numEvt(BdkEvtTypes::SIG_GOOD_D);
    RooRealVar * asym = pdfOnResDK.typeAsym(BdkEvtTypes::SIG_GOOD_D);
    
    cout << "Changing nsig and asym from " 
	 << nsig->getVal() << ", " << asym->getVal() << " to ";
    
    nsig->setVal(rand->Gaus(nsig->getVal(), nsig->getError()));
    asym->setVal(rand->Gaus(asym->getVal(), asym->getError()));
    
    cout << nsig->getVal() << ", " << asym->getVal() << endl;
  }

  if (toy) {
    if (useShape || useBoth) {
      if (nEvents<=0) {
        cout << "Generating full toy experiment" << endl;
        data = pdfOnResDK.generate();
      }
      else {
        cout << "Generating toy experiment with " <<nEvents<<" events"<< endl;
        data = pdfOnResDK.generate(nEvents);
      }
    }
    else {
      cout << "Generating one dummy event" << endl;
      data = dalitzHolderN.sigGoodD0()->generate(1); // just so we have a pointer for fit()
    }
  }    
  else {
    cout << "Will fit external RooDataSet * data" << endl;
  }

  BdkOnResNLLYields nll(pdfOnResDK, -0.04);
  nll.Print();
  ///  nll.setVerbose("V");

  // RooMinuit.cc line 833 gives error status if the NLL returns 0. So
  // move out of the minimum a little:
  dalitzHolderN.sigGoodD0Type().x()->setVal(dalitzHolderN.sigGoodD0Type().x()->getVal() + 0.0001);

  // Refix:
  pdfOnResDK.fixAll();
  dalitzHolderN.sigGoodD0Type().x()->setConstant(kFALSE);
  dalitzHolderN.sigGoodD0Type().y()->setConstant(kFALSE);
  dalitzHolderP.sigGoodD0Type().x()->setConstant(kFALSE);
  dalitzHolderP.sigGoodD0Type().y()->setConstant(kFALSE);


  // Fit:  
  RooFitResult * res1 = 0;
  RooFitResult * res2 = 0;  
  RooFitResult * res3 = 0;

  if (useShape) {
    RooFitResult * res1 = fit(pdfOnResDK, *data);
    res1->Write();
  }

  if (useYields) {
    RooFitResult * res2 = fit(pdfOnResDK, *data, &nll, kTRUE);
    res2->Write();
  }  

  if (useBoth) {
    RooFitResult * res3 = fit(pdfOnResDK, *data, &nll);
    res3->Write();
  }


  // Summarize fit results:
  if (0 != res1) {
    cout << "=============== no penalty =============" << endl;
    printFitResult(res1);
  }
  if (0 != res2) {
    cout << "=============== penalty only =============" << endl;
    printFitResult(res2);
  }
  if (0 != res3) {
    cout << "=============== nll + penalty =============" << endl;
    printFitResult(res3);
  }
}



