// $Id: plotDalitzEvtTypes.cc,v 1.1 2006/04/05 01:51:45 fwinkl Exp $
// Makes a Dalitz plot for the different event types
// Only the B- mode is plotted.
// By default only MC truth cuts are done. Supply additional cuts in addCut.

void plotDalitzAllEvtTypes(TCut addCut="") {

  dalitzEvtTypes(sigTree,cutDKBadD&&addCut);
  dalitzEvtTypes(sigTree,cutDKGoodD&&addCut);
  dalitzEvtTypes(bbTree,cutDKX&&addCut);
  dalitzEvtTypes(bbTree,cutDPiX&&addCut); 
  dalitzEvtTypes(bbTree,cutDPiBadD&&addCut);
  dalitzEvtTypes(bbTree,cutDPiGoodD&&addCut);
  dalitzEvtTypes(bbTree,cutBBGoodD&&addCut);
  dalitzEvtTypes(bbTree,cutBBBadD&&addCut);
  dalitzEvtTypes(qqTree,cutqqGoodD&&addCut);
  dalitzEvtTypes(qqTree,cutqqBadD&&addCut);
}

void plotDalitzEvtTypes(TTree *tree, TCut cut) {
  
  gStyle->SetOptStat(1111);
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.85);
  gStyle->SetTitleYOffset(1.2);
  TCanvas *can = new TCanvas("can","Dalitz plot",500,500);
  can->SetTopMargin(0.05);
  //  can->SetRightMargin(0.05);

  const int bins = 100;
    
  TH2D *h2 = new TH2D(cut.GetName(),"",bins,0,3,bins,0,3);
  TString cmd = "d0ppmupmass**2:d0pppupmass**2>>";
  cmd += cut.GetName();

  tree->Draw(cmd,cut&&"(Hdtrkchge<0)","goff");
  h2->GetXaxis()->SetTitle("M(#pi^{+}#pi^{0})^{2}");
  h2->GetYaxis()->SetTitle("M(#pi^{-}#pi^{0})^{2}");
  h2->SetMarkerStyle(21);
  h2->SetMarkerColor(kBlue);
  h2->SetMarkerSize(0.6);
  h2->Draw();
  
  can->SaveAs(TString(cut.GetName())+".root");
  can->SaveAs(TString(cut.GetName())+".png");

}

