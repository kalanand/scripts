// $Id: smearDataSet.cc,v 1.1 2006/03/16 07:08:28 fwinkl Exp $
// Smear a dataset with a resolution function

RooDataSet* smearDataSet(RooDataSet& data, RooRealVar& m12, RooRealVar& m13, BdkPdfAbsBase& m12Res, BdkPdfAbsBase& m13Res)
{
  
  if (data.get()->find(m12.GetName())==0) {
    cout << "Dataset does not contain variable "<<m12.GetName()<<endl;
    return 0;
  }
  if (data.get()->find(m13.GetName())==0) {
    cout << "Dataset does not contain variable "<<m13.GetName()<<endl;
    return 0;
  }

  // Create smeared dataset
  RooDataSet *smear = new RooDataSet(TString(data.GetName())+".smeared",
                                     TString(data.GetTitle())+" smeared",
                                     *data.get());

  // Create resolution datasets
  RooDataSet* m12ResData = m12Res.generate(data.numEntries());
  RooDataSet* m13ResData = m13Res.generate(data.numEntries());

  if (m12ResData->get()->find(m12.GetName())==0) {
    cout << "m12 resolution PDF does not depend on variable "<<m12.GetName()<<endl;
    return 0;
  }
  if (m13ResData->get()->find(m13.GetName())==0) {
    cout << "m13 resolution PDF does not depend on variable "<<m13.GetName()<<endl;
    return 0;
  }
  
  
  for (Int_t i = 0; i<data.numEntries(); i++) {
    RooArgSet vars(*data.get(i));
    
    RooRealVar *r = (RooRealVar*)vars.find(m12.GetName());
    RooRealVar *s = (RooRealVar*)m12ResData->get(i)->find(m12.GetName());    
    r->setVal(r->getVal() + s->getVal());
    
    r = (RooRealVar*)vars.find(m13.GetName());
    s = (RooRealVar*)m13ResData->get(i)->find(m13.GetName());
    r->setVal(r->getVal() + s->getVal());
    
    smear->add(vars);
  }

  delete m12ResData;
  delete m13ResData;
  
  return smear;
}
