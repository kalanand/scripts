// $Id: weightedMCHisto.cc,v 1.2 2006/05/29 22:57:15 fwinkl Exp $
//
// fill histogram with weighted MC
// works for 1D and 2D
// for 2D use var = "yvar:xvar" notation as in the TTree::Draw() command

void weightedMCHisto(TH1* h, const char* var, TCut cut, Bool_t addSignal = kFALSE)
{
  TString cuts;

  h->Sumw2();

  if (b0Tree) {
    cuts.Form("%f*(%s)",B0_WEIGHT,cut.GetTitle());
    b0Tree->Project(h->GetName(),var,cuts);
  }
  else cout <<"weightedMCHisto(): No B0 tree. Skipping."<<endl;

  if (bpTree) {
    cuts.Form("%f*(%s)",BP_WEIGHT,cut.GetTitle());
    bpTree->Project(TString("+")+h->GetName(),var,cuts);
  }
  else cout <<"weightedMCHisto(): No B+ tree. Skipping."<<endl;

  if (dpiTree) {
    cuts.Form("%f*(%s)",DPI_WEIGHT,cut.GetTitle());
    dpiTree->Project(TString("+")+h->GetName(),var,cuts);
  }
  else cout <<"weightedMCHisto(): No DPi tree. Skipping."<<endl;

  if (udsTree) {
    cuts.Form("%f*(%s)",UDS_WEIGHT,cut.GetTitle());
    udsTree->Project(TString("+")+h->GetName(),var,cuts);
  }
  else cout <<"weightedMCHisto(): No UDS tree. Skipping."<<endl;

  if (ccTree) {
    cuts.Form("%f*(%s)",CC_WEIGHT,cut.GetTitle());
    ccTree->Project(TString("+")+h->GetName(),var,cuts);
  }
  else cout <<"weightedMCHisto(): No CC tree. Skipping."<<endl;
  
  if (addSignal) {
    if (sigTree) {
      cuts.Form("%f*(%s)",SIG_WEIGHT,cut.GetTitle());
      sigTree->Project(TString("+")+h->GetName(),var,cuts);
    }
    else cout <<"weightedMCHisto(): No signal tree. Skipping."<<endl;
  }
}
