void newChunksToy(const char* resultFile, int todo=0, int bin = 0,int bbin = 0,
               const char * moreFixFile = "analysis/defaultMoreFixFile.cc",
               const char * parInputFile = "analysis/defaultParInputFile.cc",
               const char * floatFile = "analysis/defaultFloatFile.cc") {

  //reload the path
  paramsOnResDK->readFromFile(paramsPath() + "all.par");
  cout<<" read par file from "<<paramsPath()<<"all.par"<<endl;
  // Execute commands before generating events:
  cout << "reading parInputFile " << parInputFile << endl;
  gROOT->Macro(parInputFile);
                                                                                                                                   
  // Fix everything and float some:
  pdfOnResDK.fixAll();
  cout << "executing floatFile " << floatFile << endl;
  gROOT->Macro(floatFile);
 
  delete data;

  int numSigBad = 0;
  int numSigGood = 0;
  double plusEnt = 0;
  double minusEnt = 0;
  
    if (todo & SIG_G_BIT) {
        cout<< endl << "===============  signal good sample  ============="<<endl;
        RooDataSet * data0 = 0;
        int totSig = addChunk(data0, chainRanDK, cutGoodD, 1);
        delete data0;
        numSigGood = addChunk(data, pdfOnResDK.sigGoodD0(), totSig, WEIGHT_SIG, bin);
        for(int i=0; i<numSigGood; ++i) {
     	 	 // Get the original #'s of +/- events:
     	 	 RooArgSet * kchargeset = data->get(i);
     		 RooCategory * kk = (RooCategory *) kchargeset->find(Hdtrkchge->GetName());
      	 	if(1 == kk->getIndex()) {
      	 		 ++plusEnt;
      		}
      		if(-1 == kk->getIndex()) {
       			 ++minusEnt;
      		}
    	}
    } 

    if (todo & SIG_B_BIT) {
        cout<< endl << "===============  signal bad sample  ============="<<endl;
        RooDataSet * data0 = 0;
        int totSigBad = addChunk(data0, chainRanDK, cutBadD, 1);
        delete data0;
        numSigBad = addChunk(data, pdfOnResDK.sigBadD0(), totSigBad, WEIGHT_SIG, bin);
        
    }
   
    // Dpi:
    int numDpiGood = 0;
    int numDpiBad = 0;
    if(todo & DPI_G_BIT) {
       cout<< endl<<"=================   D0 pi good ============="<<endl;
       RooDataSet * data0 = 0; 	
       int totDpiGood = addChunk(data0, chainRanBchDpi, cutGoodD, 1);
       delete data0;
       numDpiGood = addChunk(data, pdfOnResDK.DpiGoodD0(), totDpiGood, WEIGHT_BCH, bbin); 
    }
    if(todo & DPI_B_BIT) {
       cout<< endl<<"=================   D0 pi bad ============="<<endl;
       RooDataSet * data0 = 0;	
       int totDpiBad  = addChunk(data0, chainRanBchDpi, cutBadD,  1);
       delete data0;
       numDpiBad = addChunk(data, pdfOnResDK.DpiBadD0(), totDpiBad, WEIGHT_BCH, bbin);
    }

    // Dpipi:
   int numDpipi = 0;
   if( todo & DPIPIBIT) {
      cout<<"=================  D0 pi pi( D(*) pi/rho )  ============="<<endl;
      RooDataSet * data0 = 0;
      TCut DpiOthCut = "((B1decmode>99&&B1decmode<109)||B1decmode>149)||((B2decmode<109&&B2decmode>99)||B2decmode>149)";
      int numDpiOther = addChunk(data0, chainRanBchComb,DpiOthCut+cutBadD, 1,  bbin);
      numDpiOther += addChunk(data0, chainRanB0Comb,DpiOthCut+cutBadD, (WEIGHT_B0/WEIGHT_BCH));
      delete data0;
      numDpipi = addChunk(data, pdfOnResDK.DPiPi(), numDpiOther, WEIGHT_BCH, bbin);
   }
 
    // charmless:
   int numCharmless = 0;
   if( todo & CHARMLESSBIT ) {
      cout<< endl<<"=================  Charmless( D* K*)  ============="<<endl;
      // DKother:
      RooDataSet * data0 = 0;
      TCut dkCut1 = "((B1decmode>9&&B1decmode<19)||(109<B1decmode&&B1decmode<149))||((B2decmode<19&&B2decmode>9)||(109<B2decmode&&B2decmode<149))";
      int numDKother = addChunk(data0, chainRanBchComb,dkCut1+cutBadD, 1 ,  bbin);
      numDKother += addChunk(data0, chainRanB0Comb,dkCut1+cutBadD, (WEIGHT_B0/WEIGHT_BCH));
      delete data0;
      numCharmless = addChunk(data, pdfOnResDK.charmless(), numDKother, WEIGHT_BCH, bbin);
   }

  // BBcomb:
  int numBBcomBad = 0;
  int numBBcomGood = 0;
  if( todo & BB_G_BIT) {
     RooDataSet * data0 = 0;
     cout<< endl<<"=================  BB good  ============="<<endl;
     int totBBcomGood = addChunk(data0, chainRanBchComb, cutGoodD, 1);
     totBBcomGood += addChunk(data0, chainRanB0Comb,  cutGoodD,(WEIGHT_B0/WEIGHT_BCH));
     delete data0;
     numBBcomGood = addChunk(data, pdfOnResDK.BBGoodD0(), totBBcomGood, WEIGHT_BCH, bbin);	
  }
                                                                                                                                             
  if( todo & BB_B_BIT) {
    RooDataSet * data0 = 0;
    cout<< endl<<"=================  BB bad   ============="<<endl;
    TCut bbo = "B1decmode<=0&&B2decmode<=0";
    int totBBcomBad = addChunk(data0, chainRanBchComb, bbo+cutBadD, 1);
    totBBcomBad += addChunk(data0, chainRanB0Comb,  bbo+cutBadD, (WEIGHT_B0/WEIGHT_BCH));
    delete data0;
    numBBcomBad = addChunk(data, pdfOnResDK.BBBadD0(), totBBcomBad, WEIGHT_BCH, bbin);
  }

  int numQQBad = 0;
  int numQQGood = 0;
  if( todo & QQ_B_BIT) {
    RooDataSet * QQset = 0;
    cout<< endl<<"=================   comtinuum bad  ============="<<endl;
    int totQQBad = addChunk(QQset, chainRanCc,  cutBadD, WEIGHT_CC);
    totQQBad += addChunk(QQset, chainRanUds, cutBadD, WEIGHT_UDS);
    delete QQset;
    numQQBad = addChunk(data, pdfOnResDK.qqBadD0(), totQQBad);
  }
  if( todo & QQ_G_BIT) {
    cout<< endl<<"=================   comtinuum good  ============="<<endl;
    RooDataSet * data0;
    int totQQGood = addChunk(data0, chainRanCc,  cutGoodD, WEIGHT_CC);
    totQQGood += addChunk(data0, chainRanUds, cutGoodD, WEIGHT_UDS);
    delete data0;
    numQQGood = addChunk(data, pdfOnResDK.qqGoodD0(), totQQGood);
  }

                                                                                      
  double SgSym = 0;
  if(0 != (plusEnt + minusEnt)) {
      SgSym = (plusEnt - minusEnt)/(plusEnt + minusEnt);
  }

  cout<<" D0K sample: Good "<< numSigGood << " Bad "<< numSigBad << endl;
  cout<<" Dpi sample: Good "<< numDpiGood <<" Bad "<< numDpiBad << endl;
  cout<<" D(*)pi/Rho sample: " << numDpipi <<endl;
  cout<<" D(*) K(*) sample : " << numCharmless <<endl;
  cout<<" Comb sample: Good "<< numBBcomGood <<" Bad "<< numBBcomBad <<endl;
  cout<<" Cont sample: Good "<< numQQGood <<" Bad "<< numQQBad <<endl;
 
  if (0 == data || data->numEntries() == 0) {
    cout << "emptry data set" << endl;
    return;
  }

  cout<<"==============set up the initial fit value =========="<<endl;
 
  setVarFix0(paramsOnResDK, "pdfOnResDK.BBBadD0NumEvts", (double)numBBcomBad);
 
  setVarFix0(paramsOnResDK, "pdfOnResDK.charmlessNumEvts", (double)numCharmless);
 
  setVarFix0(paramsOnResDK, "pdfOnResDK.BBGoodD0Frac",
         numBBcomBad ? (double)numBBcomGood/numBBcomBad : 0);
        
  setVarFix0(paramsOnResDK, "pdfOnResDK.DpiGoodD0NumEvts", numDpiGood);
  setVarFix0(paramsOnResDK, "pdfOnResDK.DpiBadD0Frac",
             numDpiGood ? (double)numDpiBad/numDpiGood : 0);
        
  setVarFix0(paramsOnResDK, "pdfOnResDK.DPiPiNumEvts", numDpipi);
 
  setVar(paramsOnResDK, "pdfOnResDK.charmlessNumEvts", numCharmless);
 
  setVarFix0(paramsOnResDK, "pdfOnResDK.qqBadD0NumEvts", numQQBad);
  setVar(paramsOnResDK, "pdfOnResDK.qqGoodD0Frac",
         numQQBad ? (double)numQQGood/numQQBad : 0);
 
  setVarFix0(paramsOnResDK, "pdfOnResDK.sigGoodD0NumEvts", numSigGood);
  setVar(paramsOnResDK, "pdfOnResDK.sigBadD0Frac",
             numSigGood ? (double)numSigBad/numSigGood : 0);
        
  setVar(paramsOnResDK, "pdfOnResDK.asymSigGoodD",SgSym);
 
  cout<<"========== end of input the intial values =========="<<endl;
 
  cout << "executing moreFixFile " << moreFixFile << endl;
  gROOT->Macro(moreFixFile);
 
  // Store initial (true) values:
  RooArgSet paramsFree = pdfOnResDK.parametersFree();
  double initValues[300];
  TIterator *iter= paramsFree.createIterator();
  RooRealVar *arg ;
  int i = -1;
  while(arg=(RooRealVar*)iter->Next()) {
    ++i;
    initValues[i] = arg->getVal();
  }
 
  // Fit:
  cout<<"Fitting the variables"<<endl;
  paramsFree.Print("V");
  RooFitResult * result = fit(pdfOnResDK, *data);
 
  // Report the results:
  ofstream ofs(resultFile);
  Int_t status = -100;
  Double_t minNll = 100;
  if(0 != result ) {
    status = result->status();
    minNll = result->minNll();
  }
  ofs << "fitStatus " << status << " " << minNll
      << " " << paramsFree.getSize();
 
 
  iter= paramsFree.createIterator();
  i = -1;
  while(arg=(RooRealVar*)iter->Next()) {
    ++i;
    double finalValue = arg->getVal();
    double error = arg->getError();
    ofs << " " << initValues[i] << " " << finalValue << " " << error << " ";
  }
  ofs << endl;
  cout << "exiting newChunksToy" << endl;
}


 
