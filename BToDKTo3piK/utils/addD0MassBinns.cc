RooAbsData * addD0MassBinns(RooAbsData * origData, 
			    RooRealVar * mD, 
			    RooCategory * D0MassReg) {
  RooAbsData * outData = new RooDataSet("outData", "outData", *origData->get());
  ((RooDataSet *) outData)->addColumn(*D0MassReg); 
  for(int i=0; i<origData->numEntries(); ++i) {
          RooArgSet * outEvt = origData->get(i);
          RooRealVar * tmpD0 = (RooRealVar *)outEvt->find(mD->GetName());
          double tmpVal = tmpD0->getVal();
          if(tmpVal<=1.825) D0MassReg->setIndex(0);
          if(tmpVal>1.825&&tmpVal<=1.850) D0MassReg->setIndex(1);
          if(tmpVal<=1.875&&tmpVal>1.850) D0MassReg->setIndex(2);
          if(tmpVal>1.875&&tmpVal<=1.900) D0MassReg->setIndex(3);
          if(tmpVal>1.900) D0MassReg->setIndex(4);
          outEvt->add(*D0MassReg);
          outData->add(*outEvt);
  }
  return outData;
} 
          


