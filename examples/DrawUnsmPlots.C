#include "TObject.h"
void DrawUnsmPlots(TString ALGO)
{
  TFile *f;
  TString ss(ALGO+"canvasOutFile.root");
  f = new TFile(ss);
  TCanvas *c;
  c = (TCanvas*)f->Get("_R80540");
  

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
  
  for(i=0;i<histList[1]->GetNbinsX();i++)
    if (histList[1]->GetBinLowEdge(i+1)>=1101)
      histList[1]->SetBinContent(i+1,0); 

  TCanvas *c = new TCanvas("c","c");
  gPad->SetLogx();
  histList[0]->SetMinimum(0.75);
  histList[0]->SetMaximum(1.02);
  histList[0]->SetTitle("");
  histList[0]->GetXaxis()->SetTitle("jet p_{T} (GeV)");
  histList[0]->GetXaxis()->SetMoreLogLabels();
  histList[0]->GetXaxis()->SetNoExponent();
  histList[0]->GetXaxis()->SetTitleFont(42);
  histList[0]->GetXaxis()->SetLabelFont(42);
  histList[0]->GetXaxis()->SetTitleSize(0.05);
  histList[0]->GetXaxis()->SetLabelSize(0.05);
  histList[0]->GetYaxis()->SetTitle("Unsmearing correction");
  histList[0]->GetYaxis()->SetTitleFont(42);
  histList[0]->GetYaxis()->SetLabelFont(42);
  histList[0]->GetYaxis()->SetTitleSize(0.05);
  histList[0]->GetYaxis()->SetLabelSize(0.05);
  histList[0]->GetYaxis()->SetNdivisions(505);

  histList[0]->SetLineWidth(2);
  histList[0]->SetLineStyle(1);
  histList[0]->SetLineColor(1);
  histList[0]->SetMarkerStyle(20);
  histList[1]->SetLineColor(2);
  histList[1]->SetMarkerColor(2);
  histList[1]->SetLineWidth(2);
  histList[1]->SetMarkerStyle(21);
  histList[1]->SetLineStyle(2);
  histList[2]->SetLineWidth(2);
  histList[2]->SetLineStyle(3);
  histList[2]->SetLineColor(4);
  histList[2]->SetMarkerColor(4);
  histList[2]->SetMarkerStyle(24);
  

  histList[0]->Draw("HIST ][");
  histList[1]->Draw("same HIST ][");
  histList[2]->Draw("same HIST ][");
  TF1 *unit = new TF1("unit","0*x+1",0,2000);
  unit->SetLineColor(1);
  unit->SetLineStyle(2);
  unit->Draw("same");

  TLegend *leg = new TLegend(0.55,0.45,0.90,0.6);
  leg->SetTextFont(42);
  leg->AddEntry(histList[0],"0.00 #leq |y| < 0.55","L");
  leg->AddEntry(histList[1],"1.10 < |y| < 1.70","L");
  leg->AddEntry(histList[2],"1.70 < |y| < 2.50","L");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  TPaveText *pave3 = new TPaveText(0.55,0.17,0.90,0.42,"NDC");
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

  c->Print(ALGO+"_Unsmearing.eps");
  
}
