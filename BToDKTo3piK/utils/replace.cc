// $Id: replace.cc,v 1.3 2006/05/29 22:57:14 fwinkl Exp $
// Some functions to replace variables in a RooDataSet

#include "RooFitCore/RooDataSet.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooSimultaneous.hh"

#include "BToDKTo3piK/BdkPdfAbsBase.hh"

#include "../BToDKTo3piK/globals/globals.hh"


// Replace vars in origData by values from newData
RooDataSet* replace(RooDataSet* origData, 
                    RooDataSet* newData,
                    const RooArgSet& vars)
{
  if (origData->numEntries() > newData->numEntries()) {
    cout << "replace: dataset "<<origData->GetName()
         << " has more entries than dataset "<<newData->GetName()<<endl;
    return 0;
  }

  // Generate the output data set:
  RooDataSet * outData = new RooDataSet("outData", "outData", *origData->get());

  TIterator* varIter = vars.createIterator();

  for (int i = 0; i < origData->numEntries(); ++i){
    // Get the original event and the toy one:
    const RooArgSet * origEvt = origData->get(i);
    const RooArgSet * newEvt = newData->get(i);
    
    // Make the new event: 
    RooArgSet outEvt(*origEvt);

    // Replace variables
    varIter->Reset();
    RooAbsArg *arg;
    while ((arg = (RooAbsArg*)varIter->Next())) {
      RooRealVar * genVar = (RooRealVar *)newEvt->find(arg->GetName());
      RooRealVar * outVar = (RooRealVar *)outEvt.find(arg->GetName());
      outVar->setVal(genVar->getVal());
    }

    // and add to the outData:
    outData->add(outEvt);
  } 

  delete varIter;

  return outData;
}





// Replaces vars in each event in the data with the pdf-generated values.
// Caller owns the return value.
RooDataSet * replace(RooDataSet * origData, 
                     const RooArgSet& vars,
		     BdkPdfAbsBase& pdfP, BdkPdfAbsBase& pdfN) {
  
  if (vars.getSize()==0 || 0 == origData) {
    return origData;
  }

  // Create sim-PDF
  RooSimultaneous sim("sim","sim",*Hdtrkchge);
  sim.addPdf(*pdfP.getPdf(), "+");
  sim.addPdf(*pdfN.getPdf(), "-");
  
  // Generate the variables to replace using the original dataset as prototype
  RooDataSet* proto = (RooDataSet*)origData->reduce(RooArgSet(*Hdtrkchge));
  RooDataSet * genData = sim.generate(vars, *proto);  

  // Check that the pdf generates all vars
  TIterator* iter = vars.createIterator();
  RooAbsArg *arg;
  while ((arg = (RooAbsArg*)iter->Next())) {    
    if (0 == genData->get(0)->find(arg->GetName())) {
      cout << "pdf does not generate "<<arg->GetName()<<". Returning null dataset." << endl;
      return 0;
    }
  }
  delete iter;
  delete proto;

  cout << "replace: replacing variables:"<<endl;
  vars.Print();
  cout << "using PDF " << pdfP.GetName() 
       << " and "<< pdfN.GetName() << endl;

  // Generate the output data set:
  RooDataSet * outData = replace(origData, genData, vars);

  delete genData;
  return outData;
}
