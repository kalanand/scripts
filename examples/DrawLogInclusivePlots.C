#include "TObject.h"
void DrawLogInclusivePlots(TString ALGO)
{
  TFile *f;
  TString ss(ALGO+"canvasOutFile.root");
  f = new TFile(ss);
  TCanvas *c
  c = (TCanvas*)f->Get("_R80570");
  TString algoTitle;
  if (ALGO=="KT6")
    algoTitle = "k_{T} D = 0.6";
  else
    algoTitle = "SISCone R = 0.7";
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
  histList[0]->SetMaximum(1e+10);
  histList[0]->SetTitle("");
  histList[0]->GetXaxis()->SetTitle("jet p_{T} (GeV)");
  histList[0]->GetXaxis()->SetTitleFont(42);
  histList[0]->GetXaxis()->SetLabelFont(42);
  histList[0]->GetXaxis()->SetTitleSize(0.05);
  histList[0]->GetXaxis()->SetLabelSize(0.05);
  histList[0]->GetXaxis()->SetMoreLogLabels();
  histList[0]->GetXaxis()->SetNoExponent();
  histList[0]->GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}dy (fb/GeV)");
  histList[0]->GetYaxis()->SetTitleFont(42);
  histList[0]->GetYaxis()->SetLabelFont(42);
  histList[0]->GetYaxis()->SetTitleSize(0.05);
  histList[0]->GetYaxis()->SetLabelSize(0.05);
  histList[0]->GetYaxis()->SetNdivisions(505);

  histList[0]->SetLineWidth(2);
  histList[0]->SetLineStyle(1);
  histList[0]->SetLineColor(2);
  histList[1]->SetLineColor(1);
  histList[1]->SetMarkerColor(1);
  histList[1]->SetLineWidth(1);
  histList[1]->SetMarkerStyle(20);
  histList[1]->SetLineStyle(1);

  histList[2]->SetLineWidth(2);
  histList[2]->SetLineStyle(1);
  histList[2]->SetLineColor(2);
  histList[3]->SetLineColor(1);
  histList[3]->SetMarkerColor(1);
  histList[3]->SetLineWidth(1);
  histList[3]->SetMarkerStyle(24);
  histList[3]->SetLineStyle(1);

  histList[4]->SetLineWidth(2);
  histList[4]->SetLineStyle(1);
  histList[4]->SetLineColor(2);
  histList[5]->SetLineColor(1);
  histList[5]->SetMarkerColor(1);
  histList[5]->SetLineWidth(1);
  histList[5]->SetMarkerStyle(21);
  histList[5]->SetLineStyle(1);

  histList[0]->Draw("C");
  histList[1]->Draw("same");
  histList[2]->Draw("sameC");
  histList[3]->Draw("same");
  histList[4]->Draw("sameC");
  histList[5]->Draw("same");

  TLegend *leg = new TLegend(0.47,0.75,0.92,0.9);
  leg->SetTextFont(42);
  leg->AddEntry(histList[1],"0.00 #leq |y| < 0.55 ( #times 32 )","P");
  leg->AddEntry(histList[3],"1.10 < |y| < 1.70 ( #times 16 )","P");
  leg->AddEntry(histList[5],"1.70 < |y| < 2.50 ( #times 8 )","P");
  leg->AddEntry(histList[0],"Theory","L");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  TPaveText *pave3 = new TPaveText(0.2,0.2,0.45,0.45,"NDC");
  pave3->SetTextFont(42);
  pave3->AddText("CMS preliminary");
  pave3->AddText(algoTitle);
  pave3->AddText("#sqrt{s} = 10 TeV");
  pave3->AddText("");
  pave3->AddText("#int L dt = 10 pb^{-1}");
  pave3->SetFillColor(0);
  pave3->SetBorderSize(0);
  pave3->SetTextSize(0.04);
  pave3->Draw();

  c->Print(ALGO+"_LogLogComparison.eps");
}
