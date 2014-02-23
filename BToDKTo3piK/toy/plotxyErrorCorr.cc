// $Id: plotxyErrorCorr.cc,v 1.2 2006/06/23 20:45:50 fwinkl Exp $
// Make scatter plots of x/y and their errors and pulls

void plotxyErrorCorr(const char* fitMC,
                     RooRealVar* xvar = dalitzHolderN.sigGoodD0Type().x(),
                     RooRealVar* yvar = dalitzHolderN.sigGoodD0Type().y())
{
  TFile* f = new TFile(fitMC);
  RooDataSet* d = (RooDataSet*)f->Get("fitParData");

  gStyle->SetOptStat(0);  
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetTitleOffset(1,"x");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.06,"x");
  gStyle->SetLabelSize(0.06,"y");
  gStyle->SetNdivisions(506,"x");
  gStyle->SetNdivisions(506,"y");

  TCanvas* can = new TCanvas("can","can",1200,325);
  can->Divide(4,1,0.001,0.001);
  can->cd(1);

  TH2D* hx = new TH2D("hx",";x;error x",100,-1,1,100,0,0.8);
  hx->Draw();
  d->tree().Draw(xvar->GetName()+TString("err:")+xvar->GetName(),"","same");

  can->cd(2);
  TH2D* hy = new TH2D("hy",";y;error y",100,-1,1,100,0,0.8);
  hy->Draw();
  d->tree().Draw(yvar->GetName()+TString("err:")+yvar->GetName(),"","same");

  can->cd(3);
  TH2D* hxy1 = new TH2D("hxy1",";x;error y",100,-1,1,100,0,0.8);
  hxy1->Draw();
  d->tree().Draw(yvar->GetName()+TString("err:")+xvar->GetName(),"","same");

  can->cd(4);
  TH2D* hxy2 = new TH2D("hxy2",";y;error x",100,-1,1,100,0,0.8);
  hxy2->Draw();
  d->tree().Draw(xvar->GetName()+TString("err:")+yvar->GetName(),"","same");

  can->SaveAs("plotxyErrorCorr.eps");  


  TCanvas* can2 = new TCanvas("can2","can2",1200,325);
  can2->Divide(4,1,0.001,0.001);

  can2->cd(1);
  TH2D* hpullx = new TH2D("hpullx",";error x;pull x",100,0,0.8,100,-5,5);
  hpullx->Draw();
  d->tree().Draw(xvar->GetName()+TString("pull:")+xvar->GetName()+TString("err"),"","same");

  can2->cd(2);
  TH2D* hpully = new TH2D("hpully",";error y;pull y",100,0,0.8,100,-5,5);
  hpully->Draw();
  d->tree().Draw(yvar->GetName()+TString("pull:")+yvar->GetName()+TString("err"),"","same");

  can2->cd(3);
  TH2D* hpullxy = new TH2D("hpullxy",";error x;pull y",100,0,0.8,100,-5,5);
  hpullxy->Draw();
  d->tree().Draw(yvar->GetName()+TString("pull:")+xvar->GetName()+TString("err"),"","same");

  can2->cd(4);
  TH2D* hpullxy2 = new TH2D("hpullxy2",";error y;pull x",100,0,0.8,100,-5,5);
  hpullxy2->Draw();
  d->tree().Draw(xvar->GetName()+TString("pull:")+yvar->GetName()+TString("err"),"","same");

  can2->SaveAs("plotPullErrorCorr.eps");  


  TCanvas* can3 = new TCanvas("can3","can3",1200,325);
  can3->Divide(4,1,0.001,0.001);

  can3->cd(1);
  TH2D* hpull2x = new TH2D("hpull2x",";x;pull x",100,-1,1,100,-5,5);
  hpull2x->Draw();
  d->tree().Draw(xvar->GetName()+TString("pull:")+xvar->GetName(),"","same");

  can3->cd(2);
  TH2D* hpull2y = new TH2D("hpull2y",";y;pull y",100,-1,1,100,-5,5);
  hpull2y->Draw();
  d->tree().Draw(yvar->GetName()+TString("pull:")+yvar->GetName(),"","same");

  can3->cd(3);
  TH2D* hpull2xy = new TH2D("hpull2xy",";x;pull y",100,-1,1,100,-5,5);
  hpull2xy->Draw();
  d->tree().Draw(yvar->GetName()+TString("pull:")+xvar->GetName(),"","same");

  can3->cd(4);
  TH2D* hpull2xy2 = new TH2D("hpull2xy2",";y;pull x",100,-1,1,100,-5,5);
  hpull2xy2->Draw();
  d->tree().Draw(xvar->GetName()+TString("pull:")+yvar->GetName(),"","same");

  can3->SaveAs("plotxyPullCorr.eps");  

  
}
