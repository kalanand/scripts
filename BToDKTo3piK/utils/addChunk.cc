// $Id: addChunk.cc,v 1.4 2006/04/18 17:57:23 fwinkl Exp $
// adds a chunk of events from tree to dat with weight

#ifndef ADDCHUNK_CC
#define ADDCHUNK_CC

#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooRealVar.hh"

#include "../BToDKTo3piK/globals/cuts.hh"
#include "../BToDKTo3piK/globals/vars.hh"
#include "../BToDKTo3piK/utils/randRep.cc"
#include "../BToDKTo3piK/utils/read.cc"

#include "BToDKTo3piK/BdkPdfAbsBase.hh"

// If dat !=0 then append data to existing dataset
int addChunk(RooDataSet *& dat, 
	     TTree* tree, 
	     TCut & cut, 
	     double weight,
	     int bin = 0) {

  readCut = cut;
  RooDataSet * origData = read(tree);
  RooDataSet * weightedData = (RooDataSet*)randRep(origData, weight, bin);
  int result = weightedData->numEntries();
  
  if (0 == dat) {
    dat = weightedData;
    return result;
  }

  ((RooDataSet*)dat)->append(*weightedData);
  return result;
}

//-------------------------------------------------
int addChunk(RooDataSet *& dat, RooDataSet * origD, 
 	    TCut cut="", double weight=1, int bin=0) {

  RooDataSet* imData = 0;
  if(strlen(cut.GetTitle())==0) {
    imData = origD;
  } else {
    imData = (RooDataSet*)origD->reduce(cut);
  }
  
  RooDataSet * weightedData = randRep(imData, weight, bin);
  int result = weightedData->numEntries();
  if (0 == dat) {
    dat = weightedData;
    return result;
  }
  
  ((RooDataSet*)dat)->append(*weightedData);
  return result;
}

//----------------------
int addChunk(RooDataSet *& dat,
             TTree* tree,
             TCut & cut,
             double weight,
             double lowcut,
             int bin = 0) {

  readCut = cut;
  RooDataSet * origData = read(tree);
  RooDataSet * weightedData = randRep(origData, (lowcut+weight), lowcut, bin);
  int result = weightedData->numEntries();

  if (0 == dat) {
    dat = weightedData;
    return result;
  }

  dat->append(*weightedData);
  return result;
}


//-----------------------

int addChunk(RooDataSet *& dat, BdkPdfAbsBase * pdf, double weight, int numEvnts) 
{
     int minEvts = int(numEvnts/2+0.5);
     RooDataSet * origData = pdf->generate(minEvts);
     Hdtrkchge->setIndex(-1);
     origData->addColumn(*Hdtrkchge);
     RooDataSet * origDataP = pdf->generate(minEvts);
     Hdtrkchge->setIndex(1); 
     origDataP->addColumn(*Hdtrkchge);
     origData->append(*origDataP);    
    // Add randRep if it doesn't already exist:
      RooDataSet * outData =
                  new RooDataSet("outData", "outData", *origData->get());
    if (0 == origData->get()->find(randAdd->GetName())) {
        outData->addColumn(*randAdd);
        TRandom ran;
        for (int i = 0; i < origData->numEntries(); ++i){
          // Get the original event:
          RooArgSet EvtSet(*origData->get(i));
          EvtSet.add(*randAdd);
          RooRealVar * tmpRan = (RooRealVar *)EvtSet.find(randAdd->GetName());
          tmpRan->setVal( ran.Rndm() );
          outData->add(EvtSet);
        }
    }
   RooDataSet * weightedData = randRep(outData, weight);
   int result = weightedData->numEntries();
   if( 0 == dat ) {
        dat = weightedData;
        return result;
   }
 
   ((RooDataSet *)dat)->append(*weightedData);
   return result;

}  	 				 

//-----------------------

int addChunk(RooDataSet *& dat, BdkPdfAbsBase * pdf, 
              int numEvnts, double weight=1,
              int bin=0 )
{
     int minEvts = int(numEvnts/2+0.5);
     RooDataSet * origData = pdf->generate(minEvts);
     Hdtrkchge->setIndex(-1);
     origData->addColumn(*Hdtrkchge);
     RooDataSet * origDataP = pdf->generate(minEvts);
     Hdtrkchge->setIndex(1);
     origDataP->addColumn(*Hdtrkchge);
     origData->append(*origDataP);
    // Add randRep if it doesn't already exist:
      RooDataSet * outData =
                  new RooDataSet("outData", "outData", *origData->get());
    if (0 == origData->get()->find(randAdd->GetName())) {
        outData->addColumn(*randAdd);
        TRandom ran;
        for (int i = 0; i < origData->numEntries(); ++i){
          // Get the original event:
          RooArgSet EvtSet(*origData->get(i));
          EvtSet.add(*randAdd);	
          RooRealVar * tmpRan = (RooRealVar *)EvtSet.find(randAdd->GetName());
          tmpRan->setVal( ran.Rndm() );
          outData->add(EvtSet);
        }
    }
   
   RooDataSet * weightedData = randRep(outData, weight, bin);
   int result = weightedData->numEntries(); 
   if( 0 == dat ) {
        dat = weightedData;
        return result;  
   } 
     
    ((RooDataSet *)dat)->append(*weightedData);
    return result;	
}

#endif
