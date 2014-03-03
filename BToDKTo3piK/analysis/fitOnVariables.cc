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

void fitOnVariables(RooRealVar * var1,int todo = 0, int bin = 0,  
                RooRealVar * var2 = 0,   int bbin=0,  
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

  TCut goodD = cutGoodD;  
  TCut badD = cutBadD; 

/*  if( mes != var1 && mes != var2 ) {
	goodD += "mes>5.265";
        badD += "mes>5.265";
  } 
*/
  if (todo & SIG_G_BIT) {
    cout<< endl << "===============  signal good sample  ============="<<endl; 
    RooDataSet * data0 = 0;
    numSigGood = addChunk(data0, chainRanDK, goodD, WEIGHT_SIG, bin);
    
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
    numSigBad = addChunk(data, chainRanDK, badD, WEIGHT_SIG, bin);
  }
  // Dpi:
  int numDpiGood = 0;
  int numDpiBad = 0;
  if(todo & DPI_G_BIT) {
    cout<< endl<<"=================   D0 pi good ============="<<endl;
    numDpiGood += addChunk(data, chainRanBchDpi, goodD, WEIGHT_BCH, bbin);
  }
  if(todo & DPI_B_BIT) {
    cout<< endl<<"=================   D0 pi bad ============="<<endl;
    numDpiBad  += addChunk(data, chainRanBchDpi, badD,  WEIGHT_BCH, bbin);
  }

  // Dpipi:
  int numDpipi = 0;
/*  if( todo & DPIPIBIT) {
    cout<<"=================  D0 pi pi  ============="<<endl;
    numDpipi += addChunk(data, chainRanBchDpipi, goodD, WEIGHT_BCH, bbin);
    numDpipi += addChunk(data, chainRanB0Dpipi,  goodD, WEIGHT_B0, bbin);
  } 
*/
   
  // charmless:
  int numCharmless = 0;
/*  if( todo & CHARMLESSBIT) {
    cout<< endl<<"=================  Charmless  ============="<<endl;
    numCharmless += addChunk(data, chainRanBchCharmless, badD, WEIGHT_BCH, bbin);
  } 
*/
// BBcomb:
  int numBBcomBad = 0;
  int numBBcomGood = 0;
  if( todo & BB_G_BIT) {
    cout<< endl<<"=================  BB good  ============="<<endl;
    numBBcomGood += addChunk(data, chainRanBchComb, goodD, WEIGHT_BCH, bbin);
    numBBcomGood += addChunk(data, chainRanB0Comb,  goodD, WEIGHT_B0, bbin);
  } 

  if( todo & BB_B_BIT) {
    cout<< endl<<"=================  BB bad   ============="<<endl;
    RooDataSet * datan = 0; 
    RooDataSet * datar = 0; 
    numBBcomBad += addChunk(data, chainRanBchComb, badD, WEIGHT_BCH, bbin);
    numBBcomBad += addChunk(data, chainRanB0Comb,  badD, WEIGHT_B0, bbin);
/*    datar = (RooDataSet *)replace(datan, mes, pdfOnResDK.BBBadD0());
    ((RooDataSet *) data)->append(*datar);
*/
    TCut dummyCut;

    numBBcomBad += addChunk(data, chainRanBchCharmless, dummyCut, WEIGHT_BCH, bbin);
    numBBcomBad += addChunk(data, chainRanBchDpipi, dummyCut, WEIGHT_BCH, bbin);
//    numBBcomBad += addChunk(data, chainRanBchDpi, badD,  WEIGHT_BCH, bbin);

  }

  int cbin = 0;

// Continuum:
  int numQQBad = 0;
  int numQQGood = 0;
  if( todo & QQ_G_BIT) {
    cout<< endl<<"=================   comtinuum good  ============="<<endl;
    numQQGood += addChunk(data, chainRanCc,  goodD, WEIGHT_CC, cbin);
    numQQGood += addChunk(data, chainRanUds, goodD, WEIGHT_UDS, cbin);
  }
  if( todo & QQ_B_BIT) {
    cout<< endl<<"=================   comtinuum bad  ============="<<endl;    
    numQQBad += addChunk(data, chainRanCc,  badD, WEIGHT_CC, cbin);
    numQQBad += addChunk(data, chainRanUds, badD, WEIGHT_UDS, cbin);
  }

   double posNum=0;
   double negNum=0;  
   RooTable * table0 = data->table(*Hdtrkchge);
   table0->Print("V");
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

  setVar(paramsOnResDK, "pdfOnResDK.sigGoodD0NumEvts", numSigGood);
  setVarFix0(paramsOnResDK, "pdfOnResDK.sigBadD0Frac", 
	     numSigGood ? (double)numSigBad/numSigGood : 0);

  setVar(paramsOnResDK, "pdfOnResDK.asymSigGoodD",SgSym);
  setVar(paramsOnResDK, "pdfOnResDK.globalAsym", gSym);	 

  cout<<"========== end of input the intial values =========="<<endl;

  // fix more vars:
  cout << "executing moreFixFile " << moreFixFile << endl;
  gROOT->Macro(moreFixFile);

  if(0 != var2 ) {
     cout<< " Two Dimesional fit on "<< var1->GetName() << " "<< var2->GetName() <<endl;    
     TwoDFit(var1, var2, data);
  } else {
     cout<< " One Dimensional fit on " << var1->GetName() <<endl; 
     OneDFit(var1, data); 
     OneDFit(var1);
  }
}



