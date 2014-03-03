#include <fstream>
void fitChunksToy(int todo=0, int b1 = 0,int b2 = 0, int tms0=0, int tms1=1,
                  const char * direct = "./",  
	          const char * moreFixFile = "analysis/defaultMoreFixFile.cc",
	          const char * parInputFile = "analysis/defaultParInputFile.cc",
	          const char * floatFile = "analysis/defaultFloatFile.cc") {


  delete data; // clean up before starting
  RooDataSet * dataSig = 0;
  RooDataSet * dataDpi = 0;
  RooDataSet * dataBch = 0;
  RooDataSet * dataB0 = 0;
  RooDataSet * dataCc = 0;
  RooDataSet * dataUds = 0;
  TCut StdCut = "";

  dataSig = randRep(chainRanDK, StdCut); 
  dataDpi = randRep(chainRanBchDpi, StdCut);
  dataBch = randRep(chainRanBchComb, StdCut);
  dataB0  = randRep(chainRanB0Comb, StdCut);
  dataCc  = randRep(chainRanCc, StdCut);
  dataUds = randRep(chainRanUds, StdCut);

  int totSigBin = int(1./WEIGHT_SIG);
  int totBBbin  = int(1./WEIGHT_BCH);

  if( b2 > totSigBin ) b2 = totSigBin;

  int bbin = 0; 
  int bin = 0;
  int tm = 0;
  for( tm=tms0; tm< tms1; ++tm) {
       for( bin=b1; bin<b2; ++bin ) {      
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
     
           delete data;
           int Seed = 100000*bbin + 1000*bin+10*tm+5;
           cout<<"input seed "<< Seed <<endl;
           setRandomGenSeed(Seed); 
     
           if (todo & SIG_G_BIT) {
              cout<< endl << "===============  signal good sample  ============="<<endl; 
              RooDataSet * data0 = 0;
              numSigGood = addChunk(data0, dataSig, cutGoodD, WEIGHT_SIG, bin);
              delete data0;
              numSigGood = addChunk(data, pdfOnResDK.sigGoodD0(), numSigGood); 
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
     	 RooDataSet * data0 = 0;
              numSigBad = addChunk(data0, dataSig, cutBadD, WEIGHT_SIG, bin);
              delete data0;
              numSigBad = addChunk(data, pdfOnResDK.sigBadD0(), numSigBad);
           }
         
         
           // Dpi:
           int numDpiGood = 0;
           int numDpiBad = 0;
           if(todo & DPI_G_BIT) {
              cout<< endl<<"=================   D0 pi good ============="<<endl;
              RooDataSet * data0 = 0;
              numDpiGood = addChunk(data0, dataDpi, cutGoodD, WEIGHT_BCH, bbin);
              delete data0;
              numDpiGood = addChunk(data, pdfOnResDK.DpiGoodD0(), numDpiGood);
           }
           if(todo & DPI_B_BIT) {
              cout<< endl<<"=================   D0 pi bad ============="<<endl;
              RooDataSet * data0 = 0;
              numDpiBad  = addChunk(data0, dataDpi, cutBadD,  WEIGHT_BCH, bbin);
              delete data0;
              numDpiBad  = addChunk(data, pdfOnResDK.DpiBadD0(), numDpiBad);
           }
         
           // Dpipi:
           int numDpipi = 0;
           if( todo & DPIPIBIT) {
             cout<<"=================  D0 pi pi( D(*) pi/rho )  ============="<<endl;
             RooDataSet * data0;
             TCut DpiOthCut = "((B1decmode>99&&B1decmode<109)||B1decmode>149)||((B2decmode<109&&B2decmode>99)||B2decmode>149)";  
             int numDpiOther = addChunk(data0, dataBch, DpiOthCut+cutBadD, WEIGHT_BCH,  bbin);
             numDpiOther += addChunk(data0, dataB0, DpiOthCut+cutBadD, WEIGHT_B0,  bbin); 
             delete data0;   
             numDpipi = addChunk(data, pdfOnResDK.DPiPi(),numDpiOther); 
           } 
         
           // charmless:
           int numCharmless = 0;
           if( todo & CHARMLESSBIT ) { 
             cout<< endl<<"=================  Charmless( D* K*)  ============="<<endl;
             // DKother:
             RooDataSet * data0 = 0;
             TCut dkCut1 = "((B1decmode>9&&B1decmode<19)||(109<B1decmode&&B1decmode<149))||((B2decmode<19&&B2decmode>9)||(109<B2decmode&&B2decmode<149))";
             int numDKother = addChunk(data0, dataBch,dkCut1+cutBadD,WEIGHT_BCH,  bbin);
             numDKother += addChunk(data0, dataB0,dkCut1+cutBadD, WEIGHT_B0,  bbin);
             delete data0; 
             numCharmless = addChunk(data, pdfOnResDK.charmless(), numDKother);
           } 
     
         // BBcomb:
           int numBBcomBad = 0;
           int numBBcomGood = 0;
           if( todo & BB_G_BIT) {
             RooDataSet * data0 = 0;
             cout<< endl<<"=================  BB good  ============="<<endl;
             numBBcomGood = addChunk(data0, dataBch, cutGoodD, WEIGHT_BCH, bbin);
             numBBcomGood += addChunk(data0, dataB0,  cutGoodD, WEIGHT_B0, bbin);
             delete data0;
             numBBcomGood = addChunk(data, pdfOnResDK.BBGoodD0(),numBBcomGood);
           } 
         
           if( todo & BB_B_BIT) {
             cout<< endl<<"=================  BB bad   ============="<<endl;
             TCut bbo = "B1decmode<=0&&B2decmode<=0";
             RooDataSet * data0 = 0;
             numBBcomBad = addChunk(data0, dataBch, bbo+cutBadD, WEIGHT_BCH, bbin);
             numBBcomBad += addChunk(data0, dataB0,  bbo+cutBadD, WEIGHT_B0,  bbin);
             delete data0;
             numBBcomBad = addChunk(data,pdfOnResDK.BBBadD0(), numBBcomBad); 
           }
         
           ++bbin;
     
            
           int cbin = 0;
         
         // Continuum:
           int numQQBad = 0;
           int numQQGood = 0;
           double QqAsym = 0;
         
           if( todo & QQ_B_BIT) {
             cout<< endl<<"=================   comtinuum bad  ============="<<endl;
             RooDataSet * QQset=0;
             numQQBad = addChunk(QQset, dataCc,  cutBadD, WEIGHT_CC, cbin);
             numQQBad += addChunk(QQset, dataUds, cutBadD, WEIGHT_UDS, cbin);
             delete QQset;
             numQQBad = addChunk(data, pdfOnResDK.qqBadD0(), numQQBad);
           }
      
          if( todo & QQ_G_BIT) {
             cout<< endl<<"=================   comtinuum good  ============="<<endl; 
             RooDataSet * data0 = 0;   
             numQQGood = addChunk(data0, dataCc,  cutGoodD, WEIGHT_CC, cbin);
             numQQGood += addChunk(data0, dataUds, cutGoodD, WEIGHT_UDS, cbin);
             delete data0;
             numQQGood = addChunk(data, pdfOnResDK.qqGoodD0(), numQQGood);
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
         /*  if(0 != (posNum - plusEnt + negNum - minusEnt )) {
         	gSym = (posNum+minusEnt - negNum-plusEnt)/(posNum -plusEnt + negNum - minusEnt);
           }
         */
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
           double totBB = numBBcomBad + numCharmless + numDpipi;
         
           cout<<"==============set up the initial fit value =========="<<endl;
         
         
           setVarFix0(paramsOnResDK, "pdfOnResDK.BBGoodD0Frac", 
         	totBB ? (double)numBBcomGood/totBB : 0);
           
           setVarFix0(paramsOnResDK, "pdfOnResDK.DpiGoodD0NumEvts", numDpiGood);
           setVarFix0(paramsOnResDK, "pdfOnResDK.DpiBadD0Frac", 
         	     numDpiGood ? (double)numDpiBad/numDpiGood : 0);
         
//           setVarFix0(paramsOnResDK, "pdfOnResDK.DPiPiNumEvts", numDpipi);
         
           setVar(paramsOnResDK, "pdfOnResDK.totBBNumEvts", totBB );

           setVarFix0(paramsOnResDK, "pdfOnResDK.DKOFrac",
             numDpipi ? (double)numCharmless/numDpipi : 0);
 
           setVarFix0(paramsOnResDK, "pdfOnResDK.DPOFrac",
             totBB ? (double)numDpipi/totBB : 0);
        
 
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
     
            TString fileName = direct;
            fileName ="fitChunksToy-";
            fileName += bin;
            fileName += "-";
            fileName += bbin;
            fileName += "-";
            fileName += tm;
            fileName += ".out";
     
         
           // Report the results:
           ofstream ofs(fileName);
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
           if( bbin >= totBBbin ) bbin = 0;
       }
   }

  cout << "exiting newChunks fit program" << endl;
}



