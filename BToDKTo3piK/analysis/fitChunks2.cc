#include <fstream>

// bits:
enum {
  SIG_G_BIT = 1,
  SIG_B_BIT = 2,
  DPI_G_BIT = 4,
  DPI_B_BIT = 8,
  CHARMLESSBIT = 16,
  DPIPIBIT = 32,
  BB_G_BIT = 64,
  BB_B_BIT = 128,
  QQ_G_BIT = 256,
  QQ_B_BIT = 1024
};

void fitChunks2(const char* resultFile, int todo=0, int bin = 0,int bbin = 0, 
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

  delete data; // clean up before starting

  // Read DK signal good & bad:
  int numSigBad = 0; 
  int numSigGood = 0;
  double plusEnt = 0; 
  double minusEnt = 0;
//  redinfine range
//  WEIGHT_SIG = 1;
/*  WEIGHT_SIG = WEIGHT_SIG/WEIGHT_BCH;
  WEIGHT_B0 = WEIGHT_B0/WEIGHT_BCH;
  WEIGHT_BCH = 1.0;  
*/ 

  if (todo & SIG_G_BIT) {
    cout<< endl << "===============  signal good sample  ============="<<endl; 
    RooDataSet * data0 = 0;
    numSigGood = addChunk(data0, chainRanDK, cutGoodD, WEIGHT_SIG, bin);
    
    data = data0;
    if(0 != data) {  
      RooTable * table = data->table(*Hdtrkchge);
      table->Print("V");
    }
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
    numSigBad = addChunk(data, chainRanDK, cutBadD, WEIGHT_SIG, bin);
  }


  // Dpi:
  int numDpiGood = 0;
  int numDpiBad = 0;
  if(todo & DPI_G_BIT) {
    cout<< endl<<"=================   D0 pi good ============="<<endl;
    numDpiGood += addChunk(data, chainRanBchDpi, cutGoodD, WEIGHT_BCH, bbin);
  }
  if(todo & DPI_B_BIT) {
    cout<< endl<<"=================   D0 pi bad ============="<<endl;
    numDpiBad  += addChunk(data, chainRanBchDpi, cutBadD,  WEIGHT_BCH, bbin);
  }

  // Dpipi:
  int numDpipi = 0;
  if( todo & DPIPIBIT) {
    cout<<"=================  D0 pi pi( D(*) pi/rho )  ============="<<endl;
    TCut DpiOthCut = "((B1decmode>99&&B1decmode<109)||B1decmode>149)||((B2decmode<109&&B2decmode>99)||B2decmode>149)";  
    int numDpiOther = addChunk(data, chainRanBchComb,DpiOthCut+cutBadD, WEIGHT_BCH,  bbin);
    numDpiOther += addChunk(data, chainRanB0Comb,DpiOthCut+cutBadD, WEIGHT_B0,  bbin);    
    numDpipi = numDpiOther; 
  } 

  // charmless:
  int numCharmless = 0;
  if( todo & CHARMLESSBIT ) { 
    cout<< endl<<"=================  Charmless( D* K*)  ============="<<endl;
    // DKother:
    TCut dkCut1 = "((B1decmode>9&&B1decmode<19)||(109<B1decmode&&B1decmode<149))||((B2decmode<19&&B2decmode>9)||(109<B2decmode&&B2decmode<149))";
    int numDKother = addChunk(data, chainRanBchComb,dkCut1+cutBadD,WEIGHT_BCH,  bbin);
    numDKother += addChunk(data, chainRanB0Comb,dkCut1+cutBadD, WEIGHT_B0,  bbin);
    numCharmless = numDKother; 
  } 
// BBcomb:
  int numBBcomBad = 0;
  int numBBcomGood = 0;
  if( todo & BB_G_BIT) {
    cout<< endl<<"=================  BB good  ============="<<endl;
    numBBcomGood += addChunk(data, chainRanBchComb, cutGoodD, WEIGHT_BCH, bbin);
    numBBcomGood += addChunk(data, chainRanB0Comb,  cutGoodD, WEIGHT_B0, bbin);
  } 

  if( todo & BB_B_BIT) {
    cout<< endl<<"=================  BB bad   ============="<<endl;
    TCut bbo = "B1decmode<=0&&B2decmode<=0";
    numBBcomBad += addChunk(data, chainRanBchComb, bbo+cutBadD, WEIGHT_BCH, bbin);
    numBBcomBad += addChunk(data, chainRanB0Comb,  bbo+cutBadD, WEIGHT_B0,  bbin);
  }

  int cbin = 0;

// Continuum:
  int numQQBad = 0;
  int numQQGood = 0;
  double QqAsym = 0;

  RooDataSet * QQset=0;

  if( todo & QQ_B_BIT) {
    cout<< endl<<"=================   comtinuum bad  ============="<<endl;
    numQQBad += addChunk(QQset, chainRanCc,  cutBadD, WEIGHT_CC, cbin);
    numQQBad += addChunk(QQset, chainRanUds, cutBadD, WEIGHT_UDS, cbin);
  }
  if( todo & QQ_G_BIT) {
    cout<< endl<<"=================   comtinuum good  ============="<<endl;    
    numQQGood += addChunk(data, chainRanCc,  cutGoodD, WEIGHT_CC, cbin);
    numQQGood += addChunk(data, chainRanUds, cutGoodD, WEIGHT_UDS, cbin);
  }

  if( 0 != QQset ) {
      if( 0!= data) { ((RooDataSet *)data)->append(*QQset); }
      else { data = (RooAbsData*)QQset; }
  

     double posQQ = 0;
     double negQQ = 0;

     for(int i=0; i< QQset->numEntries(); ++i) {
      // Get the original #'s of +/- events:
        RooArgSet * qkchSet = QQset->get(i);
        RooCategory * QCat = (RooCategory *) qkchSet->find(Hdtrkchge->GetName());
        if(1 == QCat->getIndex()) {
           ++posQQ;
       }
       if(-1 == QCat->getIndex()) {
           ++negQQ;
       }
    }
   if( (posQQ+negQQ)!=0 ) {
        QqAsym = (posQQ-negQQ)/(posQQ+negQQ);
   }

  }	  

   double posNum=0;
   double negNum=0;  
   if( 0 != data ) {
      RooTable * table0 = data->table(*Hdtrkchge);
      table0->Print("V");
   }
   for(int i=0; i< data->numEntries(); ++i) {
     // Get the original #'s of +/- events:
     RooArgSet * kchSet = data->get(i);
     RooCategory * kCat = (RooCategory *) kchSet->find(Hdtrkchge->GetName());
     if(1 == kCat->getIndex()) {
       ++posNum;
     }
     if(-1 == kCat->getIndex()) {
       ++negNum;
     }
   }

  

  double gSym = 0;
  if(0 != (posNum - plusEnt + negNum - minusEnt )) {
	gSym = (posNum+minusEnt - negNum-plusEnt)/(posNum -plusEnt + negNum - minusEnt);
  }
  double SgSym = 0;
  if(0 != (plusEnt + minusEnt)) {
      SgSym = (plusEnt - minusEnt)/(plusEnt + minusEnt); 
  }	

  cout<<"================ end reading and loading samples ============="<<endl;

  cout<<" D0K sample: Good "<< numSigGood << " Bad "<< numSigBad << endl;
  cout<<" Dpi sample: Good "<< numDpiGood <<" Bad "<< numDpiBad << endl;
  cout<<" Dpp sample: " << numDpipi <<endl;
  cout<<" charmless : " << numCharmless <<endl;
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

  setVarFix0(paramsOnResDK, "pdfOnResDK.asymQqBadD", QqAsym);
  
  setVarFix0(paramsOnResDK, "pdfOnResDK.sigGoodD0NumEvts", numSigGood);
  setVar(paramsOnResDK, "pdfOnResDK.sigBadD0Frac", 
	     numSigGood ? (double)numSigBad/numSigGood : 0);

  setVar(paramsOnResDK, "pdfOnResDK.asymSigGoodD",SgSym);
  setVarFix0(paramsOnResDK, "pdfOnResDK.globalAsym", gSym);	 

  cout<<"========== end of input the intial values =========="<<endl;

  // fix more vars:
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
  
  stdPlot();

  cout << "exiting fitChunks" << endl;
}



