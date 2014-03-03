// $Id: plotNLL_RhoPhases.cc,v 1.3 2006/07/11 21:42:03 fwinkl Exp $
// Plot the x/y NLL with changing rho- phase.
// All other rho phases are set to zero.

void plotNLL_RhoPhases() {

  TCanvas* c = new TCanvas("c","",1200,600);
  c->Divide(5,2);
 
  BdkPdfDKDalitz& pdf = dalitzHolderN.sigGoodD0Type();
  BdkDDalitzAmp& amp = *pdf.dalitzAmp();

  //  pdf.parameters().readFromFile("../BToDKTo3piK/params/cleo-noEff.par");

  // Set all amplitudes and phases to zero
  for (int i=0; i<amp.nComps(); i++) {
    amp.ampRes(i)->setVal(0);
    amp.phaseRes(i)->setVal(0);
  }

  // set CLEO parameters
  setVar(pdf.parameters(),"dalitzHolderN.sigGoodD0.pdf.dalitzAmp.Rho+_amp",1);
  setVar(pdf.parameters(),"dalitzHolderN.sigGoodD0.pdf.dalitzAmp.Rho-_amp",0.65);
  setVar(pdf.parameters(),"dalitzHolderN.sigGoodD0.pdf.dalitzAmp.Rho0_amp",0);

  RooRealVar* r = (RooRealVar*)pdf.parameters().find("dalitzHolderN.sigGoodD0.pdf.dalitzAmp.Rho-_phase");

  int i = 0;
  for (double x = 0; x<=180; x+=20) {

    r->setVal(x);
    data = pdf.generate(200);
  
    c->cd(++i);
    RooPlot* px = pdf.getPdf()->plotNLLOn(((RooRealVar*)pdf.x())->frame(-10,10),data);
    RooPlot* py = pdf.getPdf()->plotNLLOn(((RooRealVar*)pdf.y())->frame(-10,10),data);
    
    px->getAttLine()->SetLineColor(kBlue);
    py->getAttLine()->SetLineColor(kBlue);
    py->getAttLine()->SetLineStyle(kDashed);

    TString title;
    title.Form("#rho^{-} phase = %3.0f",x);
    px->SetTitle(title);
    py->SetTitle(title);

    px->GetXaxis()->SetTitle("x/y");    
    py->GetXaxis()->SetTitle("x/y");    
    px->GetYaxis()->SetTitle("");    
    py->GetYaxis()->SetTitle("");

    px->Draw();
    py->Draw("same");
 
    delete data;
    c->Update();
  }
  c->SaveAs("plotNLL_RhoPhases.eps");
  c->SaveAs("plotNLL_RhoPhases.root");
}

