void DrawDistributions(const char algorithm[100], int etabin, int ptbin)
{
  const int NETA = 6;
  const double eta_boundaries[NETA+1]={0.0,0.55,1.1,1.7,2.5,3.2,5.0};
  const double Pt[18]={5,15,25,35,45,60,75,90,120,150,200,300,400,550,750,1000,1500,5000};
  char name[1024];
 // gStyle->SetOptFit(000);

  sprintf(name,"InclusiveNoteResults_%s_AllJets.root",algorithm);
  TFile *inf = new TFile(name);
  sprintf(name,"Response_RefPt%d_Eta%d",ptbin,etabin);
  TDirectoryFile *dir = (TDirectoryFile*)inf->Get("FittedHistograms");
  TH1F *hResp = (TH1F*)dir->Get(name);
  TF1 *fit = (TF1*)hResp->GetFunction("g");
  fit->SetRange(0,2);
  fit->SetNpx(500);
  //TH1F *h = new TH1F("tmp","tmp",100,0,3);
  TF1 *h = new TF1("f","0*x+1e-6",0,3);

  TPaveText *pave = new TPaveText(0.55,0.7,0.95,0.9,"NDC");
  pave->AddText("Monte Carlo Truth");
  sprintf(name,"%1.2f < |y| < %1.2f",eta_boundaries[etabin],eta_boundaries[etabin+1]);
  pave->AddText(name);
  sprintf(name,"%d < particle jet p_{T} < %d GeV",Pt[ptbin],Pt[ptbin+1]);
  pave->AddText(name);
  if (strcmp(algorithm,"SC7")==0)
    pave->AddText("Seedless Cone R = 0.7");
  else
    pave->AddText("k_{T} D = 0.6");
  pave->SetFillColor(0);
  pave->SetLineColor(0);
  pave->SetBorderSize(0);
  pave->SetTextFont(42);

  TCanvas *c = new TCanvas("Response","Response");
  gPad->SetLogy();
  h->SetTitle("");
  h->SetMaximum(1.1*fit->GetParameter(0));
  h->GetXaxis()->SetTitle("Corrected jet response");
  h->GetXaxis()->SetNdivisions(505);
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(0.9);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetLabelOffset(0.007);
  h->GetYaxis()->SetTitle("Arbitrary units");
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(1.25);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelOffset(0.007);
  hResp->SetFillColor(kGray);
  h->Draw();
  hResp->Draw("same");
  pave->Draw();
}
