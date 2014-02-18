#include "TObject.h"
void DrawVariousPlots()
{
  TFile *f;
  TString ss("Trigger_canvasOutFile.root");
  f = new TFile(ss);
  TCanvas *c
  c = (TCanvas*)f->Get("_R38335");
  TString algoTitle;
  algoTitle = "k_{T} D = 0.6";
  TList *li = (TList*)gPad->GetListOfPrimitives();
  TList *li1;
  TObject *obj;
  
  TIter next(li);
  TH1D *histList[100], *htmp;
  TGraphErrors *gList[100];
  TF1 *fList[100];
  int i,j,N(0),N1(0);
  

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
          li1 = (TList*)histList[N]->GetListOfFunctions();
          TIter next1(li1);
          fList[N] = (TF1*)next1();
          N++;
        } 
      if (cname=="TGraphErrors")
        {
          gList[N1] = (TGraphErrors*)gPad->FindObject(obj);
          gList[N1]->SetFillColor(0);
          gList[N1]->SetFillStyle(0);
          N1++;
        }
    }
  
  TCanvas *c = new TCanvas("c","c");
  fList[0]->SetMaximum(1.5);
  fList[0]->SetMinimum(0.5);
  fList[0]->SetTitle("");
  fList[0]->SetRange(20,140);
  fList[0]->GetXaxis()->SetTitle("jet p_{T} (GeV)");
  fList[0]->GetXaxis()->SetTitleOffset(1.2);
  fList[0]->GetXaxis()->SetTitleFont(42);
  fList[0]->GetXaxis()->SetLabelFont(42);
  fList[0]->GetXaxis()->SetTitleSize(0.05);
  fList[0]->GetXaxis()->SetLabelSize(0.05);
  fList[0]->GetYaxis()->SetTitle("Relative Trigger Efficiency");
  fList[0]->GetYaxis()->SetTitleFont(42);
  fList[0]->GetYaxis()->SetLabelFont(42);
  fList[0]->GetYaxis()->SetTitleSize(0.05);
  fList[0]->GetYaxis()->SetLabelSize(0.05);
  fList[0]->GetYaxis()->SetNdivisions(505);

  fList[0]->SetLineWidth(4);
  fList[0]->SetLineStyle(1);
  fList[0]->SetLineColor(1);

  fList[1]->SetLineWidth(4);
  fList[1]->SetLineStyle(2);
  fList[1]->SetLineColor(2);

  fList[2]->SetLineWidth(4);
  fList[2]->SetLineStyle(3);
  fList[2]->SetLineColor(4);

  fList[3]->SetLineWidth(4);
  fList[3]->SetLineStyle(4);
  fList[3]->SetLineColor(6);

  gList[0]->SetLineColor(8);
  gList[1]->SetLineColor(8);
  gList[2]->SetLineColor(8);
  gList[3]->SetLineColor(8); 

  fList[0]->Draw("L");
  fList[1]->Draw("sameL");
  fList[2]->Draw("sameL");
  fList[3]->Draw("sameL");
  gList[0]->Draw("same >");
  gList[1]->Draw("same >");
  gList[2]->Draw("same >");
  gList[3]->Draw("same >");

  TF1 *unit = new TF1("unit","0*x+1",0,2000);
  unit->SetLineColor(1);
  unit->SetLineStyle(2);
  unit->Draw("same");

  TLegend *leg = new TLegend(0.6,0.65,0.9,0.9);
  leg->SetTextFont(42);
  leg->AddEntry(fList[0],"HLT_Jet30","L");
  leg->AddEntry(fList[1],"HLT_Jet50","L");
  leg->AddEntry(fList[2],"HLT_Jet80","L");
  leg->AddEntry(fList[3],"HLT_Jet110","L");
  leg->AddEntry(gList[0],"Turn-on point","L");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  TPaveText *pave3 = new TPaveText(0.2,0.65,0.45,0.9,"NDC");
  pave3->SetTextFont(42);
  pave3->AddText("CMS preliminary");
  pave3->AddText(algoTitle);
  pave3->AddText("#sqrt{s} = 10 TeV");
  pave3->AddText("|y| < 0.55");
  pave3->SetFillColor(0);
  pave3->SetBorderSize(0);
  pave3->Draw();

  c->Print("KT6_TurnOn.eps");
 
}
