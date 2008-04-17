// $Id: fitMC.cc,v 1.11 2006/07/03 17:37:03 fwinkl Exp $
// Script to fit the weighted MC sample
// For serious use, compile it. Otherwise it will crash on you randomly.
// Run with analysis/runFitMC.csh

#include "TString.h"
#include "TCut.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"

#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooCategory.hh"
#include "RooFitCore/Roo1DTable.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooFitResult.hh"

#include "RooFitModels/RooGaussian.hh"


#include "../BToDKTo3piK/analysis/fitShapedFlat.cc"
#include "../BToDKTo3piK/utils/read.cc"
#include "../BToDKTo3piK/utils/misc.cc"
#include "../BToDKTo3piK/utils/addChunk.cc"
#include "../BToDKTo3piK/utils/trueAsym.cc"
#include "../BToDKTo3piK/utils/replace.cc"
#include "../BToDKTo3piK/utils/setVar.cc"
#include "../BToDKTo3piK/utils/floatVar.cc"
#include "../BToDKTo3piK/globals/globals.hh"
#include "../BToDKTo3piK/globals/setupred.cc"

#include "BToDKTo3piK/BdkSplitData.hh"
#include "BToDKTo3piK/BdkDKDalitz.hh"
#include "BToDKTo3piK/BdkPdfDDalitzInc.hh"


// Global variables used by fitSample
// Tried to give them as parameters but apparently there is a limit in CINT (crashes)
RooDataSet* bpChunk = 0; Double_t bp_w = 0;
RooDataSet* b0Chunk = 0; Double_t b0_w = 0;
RooDataSet* udsChunk = 0; Double_t uds_w = 0;
RooDataSet* ccChunk = 0; Double_t cc_w = 0;
RooDataSet* dpiChunk = 0; Double_t dpi_w = 0;
RooDataSet* sigGoodBpChunk = 0; Double_t sigGoodBp_w = 0;
RooDataSet* sigGoodBmChunk = 0; Double_t sigGoodBm_w = 0;
RooDataSet* sigChunk = 0; Double_t sig_w = 0;

RooDataSet* sigGoodData = 0;
RooArgSet* initParams = 0;

// forward declarations
void fitSample(const char* baseName = 0,
               int todo = ALL_TYPES, int bin = 0,
               RooArgSet* replaceVars = 0, Int_t toreplace = 0);

void fitMC(Int_t todo = ALL_TYPES,
           const char* baseName="fitMC", Int_t nFits = -1,
           RooArgSet* replaceVars = 0, Int_t toreplace = 0,
           Double_t scale = 1.0, const char* initParFile = 0);


// This one needs to be called, if compiled
void compFitMC(Int_t todo = ALL_TYPES,
               const char* baseName="fitMC", Int_t nFits = -1,
               RooArgSet* replaceVars = 0, Int_t toreplace = 0,
               Double_t scale = 1.0, const char* initParFile = 0)
{
  setupred();
  fitMC(todo,baseName,nFits,replaceVars,toreplace,scale,initParFile);
}


// todo          what to fit
// baseName      e.g. awg/fitMC/fitMC-1
// nFits         number of fits (-1 = all)
// replaceVars   which variables to replace
// toreplace     which event types to replace
// scale         scale all weights by this number (e.g. 3.0)
// initParFile   optional .par file to read
void fitMC(Int_t todo, const char* baseName, Int_t nFits,
           RooArgSet* replaceVars, Int_t toreplace, Double_t scale,
           const char* initParFile)
{
  TCut cut = cutSigReg;

  //reload the PDF params
  useBothFitVars();
  readOnResDKPar();
  if (initParFile) pdfOnResDK.parameters().readFromFile(initParFile);

  // Scale nBB (important for Yields-NLL)
  pdfOnResDK.nBB()->setVal(scale*pdfOnResDK.nBB()->getVal());

  cout << "todo = "<<printBits(todo) << endl;
  cout << "toreplace = "<<printBits(toreplace) << endl;
  if (replaceVars) {
    cout << "Replacing vars:"<<endl;
    replaceVars->Print();
  }
  cout << "Floating parameters:"<<endl;
  pdfOnResDK.parametersFree().Print("v");

  // save them to par file
  TString parFile = TString(baseName)+".par";
  cout << "Saving copy of PDF parameters to "<<parFile<<endl;
  pdfOnResDK.parameters().writeToFile(parFile);

  // Set asymmetry and number of signal events from x/y
  pdfOnResDK.setNsigAsymFromXY();
  Double_t nSigEvents = pdfOnResDK.numEvt(BdkEvtTypes::SIG_GOOD_D)->getVal();
  Double_t sigAsym = pdfOnResDK.typeAsym(BdkEvtTypes::SIG_GOOD_D)->getVal();

  // Read the data
  readCut = cut;
  RooDataSet* sigData = read(sigTree);
  RooDataSet* bpData = read(bpTree);
  RooDataSet* b0Data = read(b0Tree);
  RooDataSet* udsData = read(udsTree);
  RooDataSet* ccData = read(ccTree);
  RooDataSet* dpiData = read(dpiTree);
  
  // Prepare split datasets
  BdkSplitData* sigSplit = new BdkSplitData("sigSplit","bad signal",*sigData);
  BdkSplitData* bpSplit = new BdkSplitData("bpSplit","B+",*bpData);
  BdkSplitData* b0Split = new BdkSplitData("b0Split","B0",*b0Data);
  BdkSplitData* udsSplit = new BdkSplitData("udsSplit","UDS",*udsData);
  BdkSplitData* ccSplit = new BdkSplitData("ccSplit","CC",*ccData);
  BdkSplitData* dpiSplit = new BdkSplitData("dpiSplit","DPI",*dpiData);

  // Set number of chunks for each dataset
  sigSplit->setNumChunks((Int_t)(1/SIG_WEIGHT/scale));
  bpSplit->setNumChunks((Int_t)(1/BP_WEIGHT/scale));
  b0Split->setNumChunks((Int_t)(1/B0_WEIGHT/scale));
  udsSplit->setNumChunks((Int_t)(1/UDS_WEIGHT/scale));
  ccSplit->setNumChunks((Int_t)(1/CC_WEIGHT/scale));
  dpiSplit->setNumChunks((Int_t)(1/DPI_WEIGHT/scale));

  // Generate good signal events from flat signal MC

  TTree* _sigFlatTree = sigFlatTree->CopyTree(cut+cutDKGoodD);
  RooDataSet* sigFlatBpData = read(_sigFlatTree,0,"Hdtrkchge>0",allVars);
  RooDataSet* sigFlatBmData = read(_sigFlatTree,0,"Hdtrkchge<0",allVars);

  RooDataSet* sigBpData = generateFromFlat(*sigFlatBpData,
                                           (BdkDKDalitz*)dalitzHolderP.sigGoodD0()->getPdf());
  RooDataSet* sigBmData = generateFromFlat(*sigFlatBmData,
                                           (BdkDKDalitz*)dalitzHolderN.sigGoodD0()->getPdf());
 
  Double_t effBp = (Double_t)sigBpData->numEntries()/sigFlatBpData->numEntries();
  Double_t effBm = (Double_t)sigBmData->numEntries()/sigFlatBmData->numEntries();
  
  cout << "Reweighting efficencies: B+ = "<<effBp<<endl
       << "                         B- = "<<effBm<<endl;    
 
  delete _sigFlatTree;
  delete sigFlatBpData;
  delete sigFlatBmData;
 
  BdkSplitData* sigBpSplit = new BdkSplitData("sigBpSplit","signal B+",*sigBpData);
  BdkSplitData* sigBmSplit = new BdkSplitData("sigBmSplit","signal B-",*sigBmData);
  //  sigBpSplit->setNumChunks((Int_t)(1/SIG_WEIGHT*effBp/scale));
  //  sigBmSplit->setNumChunks((Int_t)(1/SIG_WEIGHT*effBm/scale));
  sigBpSplit->setChunkSize((Int_t)(nSigEvents*(1-sigAsym)/2.0));
  sigBmSplit->setChunkSize((Int_t)(nSigEvents*(1+sigAsym)/2.0));

  // Put all in a list for convenience
  TList splitList;
  splitList.Add(bpSplit);
  splitList.Add(b0Split);
  splitList.Add(udsSplit);
  splitList.Add(ccSplit);
  splitList.Add(dpiSplit);
  splitList.Add(sigSplit);
  splitList.Add(sigBpSplit);
  splitList.Add(sigBmSplit);
  
  for (int i=0; i<splitList.GetSize(); i++) {
    BdkSplitData* s = (BdkSplitData*)splitList.At(i);
    cout << "Number of chunks for "<<s->GetTitle()<<": "<<s->numChunks()
         << " ("<<s->chunkSize()<<" events each)"<<endl;
  }
  
  // We only have one set worth of contiuum events
  if (udsSplit->numChunks()>0) {
    udsChunk = (RooDataSet*)udsSplit->nextChunk();
    uds_w = UDS_WEIGHT*scale / udsSplit->chunkFraction();
  }

  if (ccSplit->numChunks()>0) {
    ccChunk = (RooDataSet*)ccSplit->nextChunk();    
    cc_w = CC_WEIGHT*scale / ccSplit->chunkFraction();
  }

  // Calculate the max number of signal and BBbar samples we have for each B+/B0
  Int_t sigSamples = TMath::Min(sigBpSplit->numChunks(),sigBmSplit->numChunks());
  Int_t bbSamples = TMath::Min(bpSplit->numChunks(),b0Split->numChunks());
  
  cout << "Number of signal samples: "<<sigSamples<<endl;
  cout << "Number of BBbar samples:  "<<bbSamples<<endl;

  for (int b=0; b<bbSamples; b++) {
    int nSigFits = sigSamples/(bbSamples-b);    
    sigSamples -= nSigFits;
    
    bpChunk = (RooDataSet*)bpSplit->chunk(b);
    b0Chunk = (RooDataSet*)b0Split->chunk(b);
    dpiChunk = (RooDataSet*)dpiSplit->chunk(b);
    
    bp_w = BP_WEIGHT*scale / bpSplit->chunkFraction();
    b0_w = B0_WEIGHT*scale / b0Split->chunkFraction();
    dpi_w = DPI_WEIGHT*scale / dpiSplit->chunkFraction();
  
    RooDataSet* initValues = 0;
    RooDataSet* finalValues = 0;
    for (int i=0; i<nSigFits; i++) {
      if (nFits==0) return;
      nFits--;
      
      cout << "========================================================="<<endl;
      cout << "      Fitting: BB sample     " << b+1<<"/"<<bbSamples<<endl;
      cout << "               signal sample " << i+1<<"/"<<nSigFits<<endl;
      cout << "========================================================="<<endl;
    
      sigGoodBpChunk = (RooDataSet*)sigBpSplit->nextChunk();
      sigGoodBmChunk = (RooDataSet*)sigBmSplit->nextChunk();
      sigChunk = (RooDataSet*)sigSplit->nextChunk();
      
      sigGoodBp_w = 1;//SIG_WEIGHT*scale / effBp / sigBpSplit->chunkFraction();
      sigGoodBm_w = 1;//SIG_WEIGHT*scale / effBm / sigBmSplit->chunkFraction();
      sig_w = SIG_WEIGHT*scale / sigSplit->chunkFraction();
      
      // Fit
      fitSample(baseName, todo, 0, replaceVars, toreplace);
      

      if (fitResult) {
        if (!initValues) {
          initParams->setAttribAll("StoreError",kTRUE);
          TString name = "initValues";
          name += b;
          initValues = new RooDataSet(name,"initial values",*initParams);
        }
        initValues->add(*initParams);

        // Only store parameters that are floating (no NLL values)
        RooArgSet* tmp = (RooArgSet*)pdfOnResDK.fitResult()->selectCommon(*initParams);
        if (!finalValues) {          
          tmp->setAttribAll("StoreError",kTRUE); 
          TString name = "finalValues";
          name += b;
          finalValues = new RooDataSet(name,"final values",*tmp);
        }
        finalValues->add(*tmp);
      }

      delete sigGoodBpChunk;
      delete sigGoodBmChunk;
      delete sigChunk;
    }

    // Save fit results
    TFile* f;
    if (b==0) f = new TFile(TString(baseName)+".root","recreate");
    else f =  new TFile(TString(baseName)+".root","update");
    if (finalValues) {
      finalValues->Write();
      delete finalValues;
    }
    if (initValues) {
      initValues->Write();
      delete initValues;
    }
    f->Close();
    delete f;
    
    delete bpChunk;
    delete b0Chunk;
    delete dpiChunk;
  }

  
  // Cleanup
  for (int i=0; i<splitList.GetSize(); i++) {    
    delete (BdkSplitData*)splitList.At(i);
  }
  delete udsChunk;
  delete ccChunk;
  
  delete sigData;
  delete sigBpData;
  delete sigBmData;
  delete bpData;
  delete b0Data;
  delete udsData;
  delete ccData;
  delete dpiData;
}


// Fit the given sample
void fitSample(const char* baseName, int todo, int bin,
               RooArgSet* replaceVars, Int_t toreplace)
{

  useBothFitVars();

  // read par file to make sure each fit uses the same configuration
  TString parFile = TString(baseName)+".par";
  pdfOnResDK.parameters().readFromFile(parFile);


  RooArgSet* paramsOnResDK = new RooArgSet(pdfOnResDK.parameters());
   
  // Cleanup before starting
  delete data;
  data = 0;

  // Dpi:
  int numDpiGood = 0;
  int numDpiBad = 0;
  if (dpiChunk && (todo & DPI_G_BIT)) {
    cout<< endl<<"=================   DPi good ============="<<endl;
    RooDataSet* dpiGData = 0;
    numDpiGood = addChunk(dpiGData, dpiChunk, cutDPiGoodD, dpi_w, bin);
    if (toreplace & DPI_G_BIT) {
      BdkPdfDDalitzInc DpiGoodD0pdf("DpiGoodD0","DpiGoodD0",*m12,*m13);
      DpiGoodD0pdf.parameters().readFromFile("DpiGoodD0.par");
  
      dpiGData = replace(dpiGData, *replaceVars, 
                         pdfOnResDK.DpiGoodD0P(), pdfOnResDK.DpiGoodD0N());
      //DpiGoodD0pdf, DpiGoodD0pdf);
    }
    if (data!=0) ((RooDataSet*)data)->append(*dpiGData);
    else data = dpiGData;
  }
  if (dpiChunk && (todo & DPI_B_BIT)) {
    cout<< endl<<"=================   DPi bad ============="<<endl;
    numDpiBad  = addChunk(data, dpiChunk, cutDPiBadD, dpi_w, bin);
  }

  // DPiX:
  int numDPiX = 0;
  if( bpChunk && b0Chunk && (todo & DPIX_BIT)) {
    cout<<"=================  DPiX  ============="<<endl;
    RooDataSet* dpixData = 0;
    numDPiX = addChunk(dpixData, bpChunk, cutDPiX, bp_w, bin);
    numDPiX += addChunk(dpixData, b0Chunk, cutDPiX, b0_w, bin);    

    if (toreplace & DPIX_BIT) {
      dpixData = replace(dpixData, *replaceVars, pdfOnResDK.DPiXP(), pdfOnResDK.DPiXN());
    }
    
    if (data!=0) ((RooDataSet*)data)->append(*dpixData);
    else data = dpixData;
    
  }
  
   
  // DKX:
  int numDKX = 0;
  if(bpChunk && b0Chunk && (todo & DKX_BIT)) {
    cout<<"=================  DKX  ============="<<endl;
    numDKX = addChunk(data, bpChunk, cutDKX, bp_w, bin);
    numDKX += addChunk(data, b0Chunk, cutDKX, b0_w, bin);    
  }

  // BBcomb:
  int numBBcomBad = 0;
  int numBBcomGood = 0;
  if (bpChunk && b0Chunk &&( todo & BB_G_BIT)) {
    cout<< endl<<"=================  BB good  ============="<<endl;
    numBBcomGood = addChunk(data, bpChunk, cutBBGoodD, bp_w, bin);
    numBBcomGood += addChunk(data, b0Chunk, cutBBGoodD, b0_w, bin);
  }
  if(bpChunk && b0Chunk &&( todo & BB_B_BIT)) {
    cout<< endl<<"=================  BB bad   ============="<<endl;
    RooDataSet * bbData = 0;
    numBBcomBad = addChunk(bbData, bpChunk, cutBBBadD, bp_w, bin);
    numBBcomBad += addChunk(bbData, b0Chunk, cutBBBadD, b0_w, bin);

    
    //if (replacer._typeBit & BB_B_BIT) {
    //  bbData = replace(bbData, replacer._var, pdfOnResDK.BBBadD0());  
    //  }
    

    if (0 != bbData) {
      if (0 != data) {
	((RooDataSet*)data)->append(*bbData);
      }
      else {
	data = bbData;
      }
    }
  }
  
  // Continuum:
  int numQQBad = 0;
  int numQQGood = 0;
  if( udsChunk && ccChunk && (todo & QQ_G_BIT)) {
    cout<< endl<<"=================   continuum good  ============="<<endl;
    numQQGood = addChunk(data, udsChunk, cutqqGoodD, uds_w, bin);
    numQQGood += addChunk(data, ccChunk, cutqqGoodD, cc_w, bin);
  }
  
  // put continuum bad in a different data set, so we can measure its asymmetry:
  RooDataSet * qqData = 0;
  if(udsChunk && ccChunk && (todo & QQ_B_BIT)) {
    cout<< endl<<"=================   continuum bad  ============="<<endl;    
    numQQBad = addChunk(qqData, udsChunk, cutqqBadD, uds_w, bin);
    numQQBad += addChunk(qqData, ccChunk, cutqqBadD, cc_w, bin);

    if (toreplace & DPIX_BIT) {
      qqData = replace(qqData, *replaceVars, pdfOnResDK.qqBadD0P(), pdfOnResDK.qqBadD0N());
    }
  }
  
  double qqBAsym = 0;

  if (0 != qqData) {
    qqBAsym = trueAsym(qqData);
    // then add qqData to data:
    if (0 != data && 0 != qqData) {
      ((RooDataSet*)data)->append(*qqData);
    }
    else {
      data = qqData;
    }
  }

  Double_t gAsym = 0;
  if (data) {
    Roo1DTable * table = data->table(*Hdtrkchge);
    table->Print("V");
    gAsym = trueAsym(data);
  }
  cout << "Global asymmetry = "<<gAsym<<endl;

  // Read DK signal good & bad:
  int numSigGood = 0;
  Double_t sigAsym = 0;
  if (sigGoodBpChunk && sigGoodBmChunk && (todo & SIG_G_BIT)) {
    cout<< endl << "===============  signal good sample  ============="<<endl; 
    
    //numSigGood = addChunk(sigData, _sigTree, cutDKGoodD, SIG_WEIGHT, bin);    
    //if (replacer._typeBit & SIG_G_BIT) {
      //      data = replace(sigData, replacer._var, pdfOnResDK.sigGoodD0());  
    //}
    
    delete sigGoodData;
    sigGoodData = 0;

    numSigGood = addChunk(sigGoodData, sigGoodBpChunk, "", sigGoodBp_w, bin);
    numSigGood += addChunk(sigGoodData, sigGoodBmChunk, "", sigGoodBm_w, bin);

    if(0 != sigGoodData) {  
      Roo1DTable * table = sigGoodData->table(*Hdtrkchge);
      table->Print("V");
      sigAsym = trueAsym(sigGoodData)-gAsym;
      cout << "Signal Asymmetry = "<<sigAsym<<endl;

      if (0 != data) ((RooDataSet*)data)->append(*sigGoodData);
      else data = sigGoodData;        
    }
  }

  int numSigBad = 0; 
  if (sigChunk && (todo & SIG_B_BIT)) {
    cout<< endl << "===============  signal bad sample  ============="<<endl; 
    numSigBad = addChunk(data, sigChunk, cutDKBadD, sig_w, bin);
  }




  cout<<"================ end reading and loading samples ============="<<endl;

  cout<<" D0K sample: Good "<< numSigGood << " Bad "<< numSigBad << endl;
  cout<<" Dpi sample: Good "<< numDpiGood <<" Bad "<< numDpiBad << endl;
  cout<<" DPiX sample: " << numDPiX <<endl;
  cout<<" DKX sample: " << numDKX <<endl;
  cout<<" BB sample: Good "<< numBBcomGood <<" Bad "<< numBBcomBad <<endl;
  cout<<" qq sample: Good "<< numQQGood <<" Bad "<< numQQBad <<endl;
   
  if (0 == data || data->numEntries() == 0) {
    cout << "empty data set" << endl;
    return;
  }

  cout<<"============== setup the initial fit values =========="<<endl;

  Int_t totBBNumEvts = numDKX + numDPiX + numBBcomBad;

  setVarFix0or1(paramsOnResDK, "pdfOnResDK.qqGoodD0Frac",
                numQQBad ? (double)numQQGood/numQQBad : 0);

  setVarFix0or1(paramsOnResDK, "pdfOnResDK.BBGoodD0Frac", 
                totBBNumEvts ? (double)numBBcomGood/totBBNumEvts : 0);

  setVarFix0or1(paramsOnResDK, "pdfOnResDK.DpiBadD0Frac", 
         numDpiGood ? (double)numDpiBad/numDpiGood : 0);

  setVarFix0or1(paramsOnResDK, "pdfOnResDK.DKXFrac", 
                numDPiX ? (double)numDKX/numDPiX : 0);

  setVarFix0or1(paramsOnResDK, "pdfOnResDK.DPiXFrac", 
                totBBNumEvts ? (double)numDPiX/totBBNumEvts : 0);

  setVarFix0or1(paramsOnResDK, "pdfOnResDK.sigBadD0Frac", 
                numSigGood ? (double)numSigBad/numSigGood : 0);

  setVarFix0(paramsOnResDK, "pdfOnResDK.DpiGoodD0NumEvts", numDpiGood);
 
  setVar(paramsOnResDK, "pdfOnResDK.asymSigGoodD",sigAsym);
  if (0 == numSigGood){
    fixVar(paramsOnResDK, "pdfOnResDK.asymSigGoodD");
    fixVar(paramsOnResDK, "dalitzHolderN.sigGoodD0.x");
    fixVar(paramsOnResDK, "dalitzHolderN.sigGoodD0.y");
    fixVar(paramsOnResDK, "dalitzHolderP.sigGoodD0.x");
    fixVar(paramsOnResDK, "dalitzHolderP.sigGoodD0.y");
  }

  setVarFix0(paramsOnResDK, "pdfOnResDK.qqBadD0NumEvts", numQQBad);
  setVarFix0(paramsOnResDK, "pdfOnResDK.sigGoodD0NumEvts", numSigGood); 
  setVarFix0(paramsOnResDK, "pdfOnResDK.totBBNumEvts", totBBNumEvts);



  setVar(paramsOnResDK, "pdfOnResDK.asymQqBadD", qqBAsym - gAsym);
  if (0 == numQQBad){
    fixVar(paramsOnResDK, "pdfOnResDK.asymQqBadD");
  }

  setVar(paramsOnResDK, "pdfOnResDK.globalAsym", gAsym);	 
  
 
  cout<<"========== end of setting up the intial values =========="<<endl;


  // Store initial (true) values:
  //  initValues = (RooArgSet*)pdfOnResDK.parametersFree().snapshot(false);

  // Fit:
  parFile = TString(baseName)+"-init.par";
  cout << "Saving initial PDF parameters to "<<parFile<<endl;
  pdfOnResDK.parameters().writeToFile(parFile);
  
  // Save initial fit parameters for pull calculation
  if (initParams) delete initParams;
  initParams = (RooArgSet*)pdfOnResDK.parametersFree().snapshot(false);

  fitResult = fit(pdfOnResDK, *data, true);


  /*
  // Report the results:
  ostringstream os;

  Int_t status = -100;
  Double_t minNll = 100;
  if(0 != result ) {
    status = result->status();
    minNll = result->minNll();
  }
  os << "fitStatus " << status << " " << minNll
     << " " << initValues->getSize()<<endl;
                                                                                
  TIterator* iter = pdfOnResDK.parametersFree().createIterator();
  RooRealVar* arg;
  while(arg = (RooRealVar*)iter->Next()) {
    double initValue = ((RooRealVar*)initValues->find(arg->GetName()))->getVal();
    double finalValue = arg->getVal();
    double error = arg->getError();
    os << arg->GetName() << " " << initValue << " "
       << finalValue << " " << error <<endl;
  }

  cout << os.str();

  if (resultFile) {
    ofstream ofs(resultFile);
    ofs << os.str();
    ofs.close();
  }
  
  //  stdPlot(data);
  */
  
  //delete data;
  //delete paramsOnResDK;
}

