#include "TObject.h"
void DrawStatUncertainty()
{
  TFile *f;
  TString ss("ExpStat_PYTHIA_canvasOutFile.root");
  f = new TFile(ss);
  TCanvas *c
  c = (TCanvas*)f->Get("_R92777");
  TString algoTitle;
  algoTitle = "k_{T} D = 0.6";
  TList *li = (TList*)gPad->GetListOfPrimitives();
  TObject *obj;
  
  TIter next(li);
  TH1D *histList[100], *htmp;
  int i,j,N;
  

  while ((obj = (TObject*)next()))
    {
      TString cname = obj->ClassName();
      TString name = obj->GetName();
      
      cout << cname <<" "<<name<<endl;    
      if (cname=="TH1D")
        {
          
          histList[N] = (TH1D*)gPad->FindObject(obj);
          histList[N]->SetFillColor(0);
          histList[N]->SetFillStyle(0);
          N++;
        } 
    }
  
  TCanvas *c = new TCanvas("c","c");
  gPad->SetLogx();
  gPad->SetLogy();
  //histList[0]->SetMaximum(1e+10);
  histList[0]->SetTitle("");
  histList[0]->GetXaxis()->SetTitle("jet p_{T} (GeV)");
  histList[0]->GetXaxis()->SetTitleFont(42);
  histList[0]->GetXaxis()->SetLabelFont(42);
  histList[0]->GetXaxis()->SetTitleSize(0.05);
  histList[0]->GetXaxis()->SetLabelSize(0.05);
  histList[0]->GetYaxis()->SetTitle("Relative Statistical Uncertainty");
  histList[0]->GetYaxis()->SetTitleFont(42);
  histList[0]->GetYaxis()->SetLabelFont(42);
  histList[0]->GetYaxis()->SetTitleSize(0.05);
  histList[0]->GetYaxis()->SetLabelSize(0.05);
  histList[0]->GetYaxis()->SetNdivisions(505);

  histList[0]->SetLineWidth(1);
  histList[0]->SetLineStyle(1);
  histList[0]->SetLineColor(1);
  histList[0]->SetMarkerColor(1);
  histList[0]->SetMarkerStyle(20);
  histList[0]->SetMarkerSize(1.2); 

  histList[1]->SetLineColor(2);
  histList[1]->SetMarkerColor(2);
  histList[1]->SetLineWidth(1);
  histList[1]->SetMarkerStyle(25);
  histList[1]->SetLineStyle(1);
  histList[1]->SetMarkerSize(1.2); 

  histList[2]->SetLineWidth(1);
  histList[2]->SetLineStyle(1);
  histList[2]->SetLineColor(4);
  histList[2]->SetMarkerColor(4);
  histList[2]->SetMarkerStyle(22);
  histList[2]->SetMarkerSize(1.2); 

  histList[0]->Draw("P");
  histList[1]->Draw("sameP");
  histList[2]->Draw("sameP");

  TLegend *leg = new TLegend(0.47,0.2,0.92,0.35);
  leg->SetTextFont(42);
  leg->AddEntry(histList[0],"|y| < 0.55","P");
  leg->AddEntry(histList[1],"1.10 < |y| < 1.70","P");
  leg->AddEntry(histList[2],"1.70 < |y| < 2.50","P");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  TPaveText *pave3 = new TPaveText(0.2,0.65,0.45,0.9,"NDC");
  pave3->SetTextFont(42);
  pave3->AddText("CMS preliminary");
  pave3->AddText(algoTitle);
  pave3->AddText("#sqrt{s} = 10 TeV");
  pave3->AddText("");
  pave3->AddText("#int L dt = 10 pb^{-1}");
  pave3->SetFillColor(0);
  pave3->SetBorderSize(0);
  pave3->Draw();

  c->Print("KT6_StatUncertainty.eps");
  
}
