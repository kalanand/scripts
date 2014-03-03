// $Id: dalitzRes.cc,v 1.2 2006/03/22 22:46:50 fwinkl Exp $
// Fit the m12/m13 Dalitz plot resolution of reconstructed MC events


void dalitzRes(TTree *tree, Bool_t fit = kFALSE) {
  
  TCut cut = cutSigReg + cutDKGoodD;
  //TCut cut = cutSigReg + cutDKBadD;
  TString baseFile = "dalitzRes";
  TString parFile = baseFile+".par";
  
   // Unsquared true Dalitz variables
  RooRealVar mass12mc("d0pppmcmass","mass12mc",1);
  RooRealVar mass13mc("d0ppmmcmass","mass13mc",1);

  RooArgSet vars(*mass12, *mass13, mass12mc, mass13mc);
  RooDataSet *data = read(tree, 0, cut, &vars);

  // Resolution of squared Dalitz variables
  RooFormulaVar m12resF("m12res","","@0-(@1*@1)",RooArgList(*m12,mass12mc));
  RooFormulaVar m13resF("m13res","","@0-(@1*@1)",RooArgList(*m13,mass13mc));
  RooRealVar m12res("m12res",TString(m12->GetTitle())+" resolution",
                    -0.05,0.05,m12->getUnit());
  //-0.1,0.1,m12->getUnit());  
  RooRealVar m13res("m13res",TString(m13->GetTitle())+" resolution",
                    -0.05,0.05,m13->getUnit());
  //  -0.1,0.1,m13->getUnit());
  
  data->addColumn(m12resF);  
  data->addColumn(m13resF);

  dalitzRes12 = BdkPdfDalitzRes("dalitzRes12","m12 Dalitz resolution",m12res);
  dalitzRes13 = BdkPdfDalitzRes("dalitzRes13","m13 Dalitz resolution",m13res);
  dalitzRes12->parameters().readFromFile("../BToDKTo3piK/params/dalitzRes.par");
  dalitzRes13->parameters().readFromFile("../BToDKTo3piK/params/dalitzRes.par");

  //dalitzRes12->parameters().readFromFile(parFile);
  //dalitzRes13->parameters().readFromFile(parFile);
  
  if (fit) {
    fitOption = "mr";
    fit(*dalitzRes12,*data);
    fit(*dalitzRes13,*data);

    ofstream of;
    of.open(parFile, ios_base::out);
    dalitzRes12->parameters().writeToStream(of,false);
    dalitzRes13->parameters().writeToStream(of,false);
    of.close();
  }

  // Plotting
  TCanvas *can = new TCanvas("can","Dalitz resolution",900,400);
  gStyle->SetTitleYOffset(1.3);
  can->Divide(2,1);  
  TGaxis::SetMaxDigits(3);

  can->cd(1);
  RooPlot *pm12 = m12res.frame();
  pm12->SetTitle("");
  data->plotOn(pm12);
  dalitzRes12->getPdf()->plotOn(pm12);
  //  dalitzRes12->getPdf()->plotOn(pm12,Components("dalitzRes12_lin.pdf"),LineStyle(kDashed));
  
  pm12->Draw();

  can->cd(2);
  RooPlot *pm13 = m13res.frame();
  pm13->SetTitle("");
  data->plotOn(pm13);
  dalitzRes13->getPdf()->plotOn(pm13);
  pm13->Draw();
 
  can->SaveAs(baseFile+".root");
  can->SaveAs(baseFile+".eps");
}



