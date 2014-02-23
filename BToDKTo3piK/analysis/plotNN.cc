// $Id: plotNN.cc,v 1.4 2006/06/28 01:59:55 fwinkl Exp $
// Plot the NN variabels for signal and background

void plotNN(TTree* sig = sigTree,
	    TTree* bb = bbTree,
	    TTree* qq = qqTree)
{
  TCut cut = cutDeltaE+cutmES+cutMD+cutDtoKpi+cutKsVeto;
  Double_t dCut = -0.25;
  
  gROOT->cd();
  TTree* _sig = sig->CopyTree(cut);
  TTree* _bb = bb->CopyTree(cut);
  TTree* _qq = qq->CopyTree(cut);

  TH1D* hqsig = new TH1D("hqsig","",20,0,1);
  TH1D* hqbb = new TH1D("hqbb","",20,0,1);
  TH1D* hqqq = new TH1D("hqqq","",20,0,1);
  hqsig->GetXaxis()->SetTitle("q");
  hqqq->SetFillStyle(3003);
  hqqq->SetFillColor(kBlue);

  _sig->Project("hqsig","nnout");
  _bb->Project("hqbb","nnout");
  _qq->Project("hqqq","nnout");

  TH1D* hdsig = new TH1D("hdsig","",20,0,1);
  TH1D* hdbb = new TH1D("hdbb","",20,0,1);
  TH1D* hdqq = new TH1D("hdqq","",20,0,1);
  
  hdsig->GetXaxis()->SetTitle("d");
  hdbb->SetFillStyle(3003);
  hdbb->SetFillColor(kBlue);

  _sig->Project("hdsig","bknnout");
  _bb->Project("hdbb","bknnout");
  _qq->Project("hdqq","bknnout");

  TCanvas* can = new TCanvas("can","Neural Network",1200,300);
  gStyle->SetTitleOffset(1.0,"x");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.06,"x");
  gStyle->SetLabelSize(0.06,"y");
  gStyle->SetNdivisions(508,"x");
  gStyle->SetNdivisions(506,"y");
  gStyle->SetOptStat(0);

  can->Divide(4,1,0.002,0.002);
  can->cd(1);
  hqsig->DrawNormalized();
  //  hqbb->DrawNormalized("same");
  hqqq->DrawNormalized("same");

  can->cd(2);
  hdsig->DrawNormalized();
  hdbb->DrawNormalized("same");
  //  hdqq->DrawNormalized("same");
  TArrow a;
  if (dCut>0) {
    a.SetFillColor(kBlack);
    a.DrawArrow(dCut,0.14,dCut,0.07,0.02);
  }

  // Plot cut efficiencies
  RooRealVar mynnout(*nnout);
  RooRealVar mybknnout(*bknnout);
  mynnout.setRange(0,1);
  mynnout.setBins(10);
  mybknnout.setRange(0,1);
  mybknnout.setBins(10);

  can->cd(3);
  plotNNeff(_sig,_qq,mynnout);

  can->cd(4);
  plotNNeff(_sig,_bb,mybknnout);

  can->SaveAs("plotNN.eps");
  can->SaveAs("plotNN.root");
}

void plotNNeff(TTree* tree1 = sigTree,
               TTree* tree2 = bbTree,
               RooRealVar& var, TCut cut = "")
{
  TGraph* g1 = 0;
  TGraph* g2 = 0;

  g1 = getCutEff(tree1,var,cut);
  if (tree2) g2 = getCutEff(tree2,var,cut);

  g1->SetMarkerStyle(8);
  g1->SetMarkerSize(0.8);
  g1->Draw("apc");
  g1->GetXaxis()->SetTitle("cut "+TString(var.GetTitle())+" > x");
  g1->GetYaxis()->SetTitle("Efficiency");


  if (g2) {
    g2->SetLineColor(kBlue);
    g2->SetMarkerStyle(22);
    g2->SetMarkerSize(0.8);
    g2->SetMarkerColor(kBlue);
    g2->Draw("pc same");
  }

  /*
  TGraph* gratio = new TGraph(g1->GetN());
  for (int i=0; i<g1->GetN(); i++) {
    Double_t x, y1, y2;
    g1->GetPoint(i,x,y1);
    g2->GetPoint(i,x,y2);
    if (y1>0) gratio->SetPoint(i,x,y2/y1);
    else if (y1==0 && y2==0) gratio->SetPoint(i,x,0);
  }
  gratio->Draw("pc same");
  */
}



TGraph* getNNeff(TTree* xtree, TTree* ytree, RooRealVar& var,
                 TCut cut = "")
{
  RooDataSet* d1 = read(xtree,0,cut,new RooArgSet(var));
  RooDataSet* d2 = read(ytree,0,cut,new RooArgSet(var));

  Double_t step = (var.getMax()-var.getMin())/var.getBins();
  TGraph* g = new TGraph(var.getBins()+1);

  int i = 0;
  Double_t nd1 = d1->sumEntries();
  Double_t nd2 = d2->sumEntries();
  for (Double_t r=var.getMin(); r<=var.getMax(); r+=step) {
    TString s;
    s.Form("%s>%f",var.GetName(),r);
    g->SetPoint(i++,d1->sumEntries(s)/nd1,d2->sumEntries(s)/nd2);
  }
  
  delete d1;
  delete d2;
  return g;
}


TGraph* getCutEff(TTree* tree, RooRealVar& var, TCut cut = "")
{
  RooDataSet* d = read(tree,0,cut,new RooArgSet(var));

  Double_t step = (var.getMax()-var.getMin())/var.getBins();
  TGraph* g = new TGraph(var.getBins()+1);
  g->SetTitle("");
  g->GetXaxis()->SetTitle(var.GetTitle());
  
  int i = 0;
  Double_t nd = d->sumEntries();
  for (Double_t r=var.getMin(); r<=var.getMax(); r+=step) {
    TString s;
    s.Form("%s>%f",var.GetName(),r);
    g->SetPoint(i++,r,d->sumEntries(s)/nd);
  }
  
  delete d;
  return g;
}
