// $Id: makeDalitzDataSet.cc,v 1.1 2006/06/14 17:57:07 fwinkl Exp $
//
// Add B- and flipped B+ Dalitz variables together into one dataset

RooDataSet* makeDalitzDataSet(const RooDataSet& data)
{  
  if (data.numEntries()==0) return 0;
  
  RooArgSet vars(*data.get());
  // Remove squared variables (we recalculate them at the end)
  vars.remove(RooArgSet(*m12,*m13),false,true);

  gROOT->cd();
  RooDataSet *dataBm = data.reduce(vars,"Hdtrkchge<0");
  RooDataSet *dataBp = data.reduce(vars,"Hdtrkchge>0");
  
  for (int i=0; i<dataBp->numEntries(); i++) {
    
    RooArgSet *args = dataBp->get(i);
    RooRealVar *arg12 = (RooRealVar*)args->find(mass12->GetName());
    RooRealVar *arg13 = (RooRealVar*)args->find(mass13->GetName());
    if (arg12==0 || arg13==0) return 0;
    
    Double_t temp12 = arg12->getVal();
    arg12->setVal(arg13->getVal());
    arg13->setVal(temp12);

    dataBm->add(*args);
  }
  delete dataBp;

  dataBm->SetName(data.GetName());
  dataBm->SetTitle(data.GetTitle());
  dataBm->addColumn(*s12);
  dataBm->addColumn(*s13);
  
  return dataBm;
}
