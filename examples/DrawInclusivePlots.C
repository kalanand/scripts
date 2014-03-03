#include "TObject.h"
void DrawInclusivePlots(TString ALGO, int Ybin)
{
  TFile *f;
  TString ss(ALGO+"_trigCombined_canvasOutFile.root");
  f = new TFile(ss);
  TCanvas *c;
  TString yTitle;
  double dY,Ymax;
  if (Ybin==0)
    {
      c = (TCanvas*)f->Get("emptyHolder_R49676");//0<|y|<0.55
      yTitle = "|y| < 0.55";
      Ymax = 2.5;
      dY = 2*0.55;
    }
  else
    {
      c = (TCanvas*)f->Get("emptyHolder_R10012");//1.70<|y|<2.50
      yTitle = "1.70 < |y| < 2.50";
      Ymax = 4;
      dY = 2*(2.5-1.7);
    }
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
  char canname[1024];
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
  /////////////////// Data ///////////////////////////
  TH1D *hData = (TH1D*)histList[15]->Clone("Data");
  
  /////////////////// Upper bounds ///////////////////
  TH1D *hJESD = (TH1D*)histList[5]->Clone("hJESD");
  hJESD->Add(histList[15],-1);

  TH1D *hPtResU = (TH1D*)histList[2]->Clone("hPtResU");
  hPtResU->Add(histList[5],-1);
  TH1D *hYResU = (TH1D*)histList[1]->Clone("hYResU");
  hYResU->Add(histList[2],-1);  

  TH1D *hPDFU = (TH1D*)histList[11]->Clone("hPDFU");
  TH1D *hNPU = (TH1D*)histList[9]->Clone("hNPU");
  hNPU->Add(hPDFU,-1);  
  TH1D *hNLOU = (TH1D*)histList[6]->Clone("hNLOU");
  hNLOU->Add(hPDFU,-1); 
  hNLOU->Add(hNPU,-1);

  /////////////////// Low bounds ///////////////////
  TH1D *hJESU = (TH1D*)histList[4]->Clone("hJESU");
  hJESU->Add(histList[15],-1);
  
  TH1D *hPtResD = (TH1D*)histList[3]->Clone("hPtResD");
  hPtResD->Add(histList[4],-1);
  TH1D *hYResD = (TH1D*)histList[0]->Clone("hYResD");
  hYResD->Add(histList[3],-1);  
  hYResD->Scale(100);  
 
  TH1D *hPDFD = (TH1D*)histList[12]->Clone("hPDFD");
  TH1D *hNPD = (TH1D*)histList[8]->Clone("hNPD");
  hNPD->Add(hPDFD,-1);  
  TH1D *hNLOD = (TH1D*)histList[7]->Clone("hNLOD");
  hNLOD->Add(hPDFD,-1); 
  hNLOD->Add(hNPD,-1);
  
  /////////////////////////////////////////////////
  TH1D *hLumiD = (TH1D*)hJESD->Clone("LumiD");
  TH1D *hLumiU = (TH1D*)hJESD->Clone("LumiU");
  for(i=0;i<hLumiD->GetNbinsX();i++)
    {
      hLumiD->SetBinContent(i+1,-0.10);
      hLumiD->SetBinError(i+1,0);
      hLumiU->SetBinContent(i+1,0.10);
      hLumiU->SetBinError(i+1,0);
    }

  TH1D *hTotExpU = (TH1D*)hJESD->Clone("TotExpU");
  AddQuadratic(hTotExpU,hPtResU);
  AddQuadratic(hTotExpU,hLumiU);
  
  TH1D *hTotExpD = (TH1D*)hJESU->Clone("TotExpD");
  AddQuadratic(hTotExpD,hPtResD);
  AddQuadratic(hTotExpD,hLumiD);
  hTotExpD->Scale(-1);

  TH1D *hTotTheU = (TH1D*)hPDFU->Clone("TotTheU");
  AddQuadratic(hTotTheU,hNPU);
  AddQuadratic(hTotTheU,hNLOU);
  
  TH1D *hTotTheD = (TH1D*)hPDFD->Clone("TotTheD");
  AddQuadratic(hTotTheD,hNPD);
  AddQuadratic(hTotTheD,hNLOD);
  hTotTheD->Scale(-1);

  for(i=0;i<hTotExpD->GetNbinsX();i++)
    {
      if (hTotExpD->GetBinCenter(i+1)<100 || hTotExpD->GetBinCenter(i+1)>1400)
        {
          hTotExpD->SetBinContent(i+1,0);
          hTotExpU->SetBinContent(i+1,0);
        }
    }  
  //////////////////////////// Experimental Error /////////////////////////////////////
  TCanvas *cExp = new TCanvas("ExperimentalError","ExperimentalError");
  gPad->SetLogx();
  TPaveText *pave = new TPaveText(0.2,0.65,0.5,0.9,"NDC");
  pave->SetTextFont(42);
  pave->AddText("CMS preliminary");
  pave->AddText(algoTitle);
  pave->AddText(yTitle);
  pave->AddText("#sqrt{s} = 10 TeV");
  pave->SetFillColor(0);
  pave->SetBorderSize(0);
  pave->SetTextSize(0.04);
  
  TLegend *leg = new TLegend(0.51,0.65,0.94,0.9);
  leg->SetTextFont(42);
  leg->AddEntry(hTotExpU,"Total","F");
  leg->AddEntry(hJESD,"Jet Energy Scale","L");
  leg->AddEntry(hLumiU,"Luminosity","L");
  leg->AddEntry(hPtResU,"Jet Energy Resolution","L");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  ///////// Style //////////////////
  hTotExpD->SetLineStyle(1);
  hTotExpD->SetLineColor(5);
  hTotExpD->SetLineWidth(1);
  hTotExpD->SetFillColor(5);
  hTotExpD->SetFillStyle(1001);
  hTotExpU->SetLineStyle(1);
  hTotExpU->SetLineColor(5);
  hTotExpU->SetLineWidth(1);
  hTotExpU->SetFillColor(5);
  hTotExpU->SetFillStyle(1001);
  hJESU->SetLineStyle(2);
  hJESU->SetLineColor(1);
  hJESU->SetLineWidth(2);
  hJESD->SetLineStyle(2);
  hJESD->SetLineColor(1);
  hJESD->SetLineWidth(2);
  hPtResU->SetLineStyle(1);
  hPtResU->SetLineColor(2);
  hPtResU->SetLineWidth(2);
  hPtResU->SetFillColor(0);
  hPtResU->SetFillStyle(1001);
  hPtResD->SetLineStyle(1);
  hPtResD->SetLineColor(2);
  hPtResD->SetLineWidth(2);
  hPtResD->SetFillColor(0);
  hPtResD->SetFillStyle(1001);
  hLumiU->SetLineStyle(3);
  hLumiU->SetLineColor(4);
  hLumiU->SetLineWidth(2);
  hLumiU->SetFillColor(0);
  hLumiD->SetLineStyle(3);
  hLumiD->SetLineColor(4);
  hLumiD->SetLineWidth(2);
  hLumiD->SetFillColor(0);
  //////////////////////////////////
  hTotExpU->SetMaximum(Ymax);
  hTotExpU->SetMinimum(-1);
  hTotExpU->SetTitle("");
  hTotExpU->GetXaxis()->SetMoreLogLabels();
  hTotExpU->GetXaxis()->SetNoExponent();
  hTotExpU->GetXaxis()->SetTitle("jet p_{T} (GeV)");
  hTotExpU->GetXaxis()->SetLabelOffset(0.015);
  hTotExpU->GetXaxis()->SetTitleFont(42);
  hTotExpU->GetXaxis()->SetLabelFont(42);
  hTotExpU->GetXaxis()->SetTitleSize(0.05);
  hTotExpU->GetXaxis()->SetLabelSize(0.05);
  hTotExpU->GetYaxis()->SetTitle("Fractional Uncertainty");
  hTotExpU->GetYaxis()->SetTitleFont(42);
  hTotExpU->GetYaxis()->SetLabelFont(42);
  hTotExpU->GetYaxis()->SetTitleSize(0.05);
  hTotExpU->GetYaxis()->SetLabelSize(0.05);
  hTotExpU->GetYaxis()->SetNdivisions(505);
  hTotExpU->Draw("AH HIST");
  hTotExpD->Draw("same HIST");
  hJESU->Draw("same HIST ][");
  hJESD->Draw("same HIST ][");
  hPtResU->Draw("same HIST ]["); 
  hPtResD->Draw("same HIST ][");
  histList[14]->Draw("same ]["); 
  hLumiU->Draw("same HIST ][");
  hLumiD->Draw("same HIST ][");
  hTotExpU->Draw("AXIS same HIST");
  leg->Draw();
  pave->Draw(); 
  sprintf(canname,"ExpError_Y%d.eps",Ybin);
  cExp->Print(ALGO+"_"+canname);  

  //////////////////////////// Theory Error /////////////////////////////////////
  TCanvas *cTheory = new TCanvas("TheoryError","TheoryError");
  gPad->SetLogx();
  TPaveText *pave1 = new TPaveText(0.2,0.65,0.5,0.9,"NDC");
  pave1->SetTextFont(42);
  pave1->AddText("CMS preliminary");
  pave1->AddText(algoTitle);
  pave1->AddText(yTitle);
  pave1->AddText("#sqrt{s} = 10 TeV");
  pave1->SetFillColor(0);
  pave1->SetBorderSize(0);
  pave1->SetTextSize(0.04);
  
  TLegend *leg = new TLegend(0.51,0.65,0.94,0.9);
  leg->SetTextFont(42);
  leg->AddEntry(hTotTheU,"Total","F");
  leg->AddEntry(hPDFU,"PDF CTEQ6.5","L");
  leg->AddEntry(hNPU,"Non-Pert. Corrections","L");
  leg->AddEntry(hNLOU,"Scale (NLO)","L");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  ///////// Style //////////////////
  hTotTheD->SetLineStyle(1);
  hTotTheD->SetLineColor(5);
  hTotTheD->SetLineWidth(1);
  hTotTheD->SetFillColor(5);
  hTotTheD->SetFillStyle(1001);
  hTotTheU->SetLineStyle(1);
  hTotTheU->SetLineColor(5);
  hTotTheU->SetLineWidth(1);
  hTotTheU->SetFillColor(5);
  hTotTheU->SetFillStyle(1001);
  hPDFU->SetLineStyle(2);
  hPDFU->SetLineColor(1);
  hPDFU->SetLineWidth(2);
  hPDFD->SetLineStyle(2);
  hPDFD->SetLineColor(1);
  hPDFD->SetLineWidth(2);
  hNLOU->SetLineStyle(3);
  hNLOU->SetLineColor(4);
  hNLOU->SetLineWidth(2);
  hNLOU->SetFillColor(0);
  hNLOU->SetFillStyle(1001);
  hNLOD->SetLineStyle(3);
  hNLOD->SetLineColor(4);
  hNLOD->SetLineWidth(2);
  hNLOD->SetFillColor(0);
  hNLOD->SetFillStyle(1001);
  hNPU->SetLineStyle(1);
  hNPU->SetLineColor(2);
  hNPU->SetLineWidth(2);
  hNPU->SetFillColor(0);
  hNPU->SetFillStyle(1001);
  hNPD->SetLineStyle(1);
  hNPD->SetLineColor(2);
  hNPD->SetLineWidth(2);
  hNPD->SetFillColor(0);
  hNPD->SetFillStyle(1001);
  //////////////////////////////////
  hTotTheU->SetMaximum(0.5);
  hTotTheU->SetMinimum(-0.2);
  hTotTheU->SetTitle("");
  hTotTheU->GetXaxis()->SetMoreLogLabels();
  hTotTheU->GetXaxis()->SetNoExponent();
  hTotTheU->GetXaxis()->SetTitle("jet p_{T} (GeV)");
  hTotTheU->GetXaxis()->SetTitleFont(42);
  hTotTheU->GetXaxis()->SetLabelFont(42);
  hTotTheU->GetXaxis()->SetLabelOffset(0.015);
  hTotTheU->GetXaxis()->SetTitleSize(0.05);
  hTotTheU->GetXaxis()->SetLabelSize(0.05);
  hTotTheU->GetYaxis()->SetTitle("Fractional Uncertainty");
  hTotTheU->GetYaxis()->SetTitleFont(42);
  hTotTheU->GetYaxis()->SetLabelFont(42);
  hTotTheU->GetYaxis()->SetTitleSize(0.05);
  hTotTheU->GetYaxis()->SetLabelSize(0.05);
  hTotTheU->GetYaxis()->SetNdivisions(505);
  hTotTheU->Draw("AH");
  hTotTheD->Draw("same");
  hNPU->Draw("same ][");
  hNPD->Draw("same ][");
  hPDFU->Draw("same ][");
  hPDFD->Draw("same ][");
  hNLOU->Draw("same ]["); 
  hNLOD->Draw("same ][");
  histList[14]->Draw("same ]["); 
  hTotTheU->Draw("AXIS same");
  leg->Draw();
  pave1->Draw();
  sprintf(canname,"TheoryError_Y%d.eps",Ybin);
  cTheory->Print(ALGO+"_"+canname);

  ////////////////////// Resolution Uncertainty //////////////
  double resErrorPlus[14] = {8.6,7.6,6.6,5.4,5.5,5.4,5.3,4.5,4.5,4.8,3.9,4.1,4.6,5.3};
  double resErrorMinus[14] = {-4.6,-4.6,-4.1,-4.5,-4.5,-4.3,-4.2,-4.5,-4.5,-4.8,-3.9,-4.1,-4.6,-5.3};
  double resErrorPlusE[6] = {6.5,5.8,5.4,5.0,5.4,4.1};
  double resErrorMinusE[6] = {-4.1,-3.8,-4.3,-3.8,-4.1,-4.1};
  double ptB[14]  = {90,112,142,175,202,230,260,292,325,357,412,525,700,1000};
  double ptE[6]  = {102,157,217,277,350,495};
  TGraph *errU = new TGraph(14,ptB,resErrorPlus);
  TGraph *errD = new TGraph(14,ptB,resErrorMinus);
  TGraph *errUE = new TGraph(6,ptE,resErrorPlusE);
  TGraph *errDE = new TGraph(6,ptE,resErrorMinusE);
  TF1 *errPlus = new TF1("errPlus","0*x+10",0,2000);
  TF1 *errMinus = new TF1("errMinus","0*x-10",0,2000);
  TCanvas *cres = new TCanvas("cres","cres");
  gPad->SetGridy();
  errU->SetMaximum(25);
  errU->SetMinimum(-15); 
  errU->SetTitle("");
  errU->GetXaxis()->SetTitle("particle jet p_{T} (GeV)");
  errU->GetYaxis()->SetTitle("Fractional Systematic Uncertainty [%]");
  errU->GetYaxis()->SetNdivisions(505);
  errU->SetLineStyle(2);
  errD->SetLineStyle(2);
  errU->SetMarkerStyle(24);
  errD->SetMarkerStyle(24);
  errU->SetMarkerSize(1.5);
  errD->SetMarkerSize(1.5);
  errUE->SetLineStyle(1);
  errDE->SetLineStyle(1);
  errUE->SetMarkerStyle(21);
  errDE->SetMarkerStyle(21);
  errUE->SetMarkerSize(1.5);
  errDE->SetMarkerSize(1.5);
  errPlus->SetLineColor(1);
  errPlus->SetLineWidth(2);
  errMinus->SetLineColor(1);
  errMinus->SetLineWidth(2);
 
  errU->Draw("ALP");
  errD->Draw("LPsame");
  errUE->Draw("LPsame");
  errDE->Draw("LPsame");
  errPlus->Draw("same");
  errMinus->Draw("same"); 

  TLegend *leg = new TLegend(0.62,0.75,0.94,0.9);
  leg->AddEntry(errU,"|#eta| < 1.1","LP");
  leg->AddEntry(errUE,"1.6 < |#eta| < 2.7","LP");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->Draw();

  TPaveText *pave2 = new TPaveText(0.2,0.75,0.6,0.9,"NDC");
  pave2->SetTextFont(42);
  pave2->AddText("Jet Energy Resolution");
  pave2->AddText("SISCone R = 0.7");
  
  pave2->SetFillColor(0);
  pave2->SetBorderSize(0);
  pave2->Draw();

  ////////////////////// (Data - Theory)/Theory //////////////
  TCanvas *cdata = new TCanvas("cdata","cdata");
  gPad->SetLogx();
  TH1D *hFinalTheoryU = (TH1D*)hTotTheU->Clone("FinalTheoryU");
  TH1D *hFinalTheoryD = (TH1D*)hTotTheD->Clone("FinalTheoryD");
  TH1D *hFinalExpU = (TH1D*)hTotExpU->Clone("FinalExpU");
  TH1D *hFinalExpD = (TH1D*)hTotExpD->Clone("FinalExpD");
  hFinalExpU->SetMaximum(2.5);
  hFinalExpU->SetMinimum(-1.5);
  hFinalExpU->GetXaxis()->SetMoreLogLabels();
  hFinalExpU->GetXaxis()->SetNoExponent();
  hFinalExpU->GetXaxis()->SetLabelOffset(0.015);
  hFinalExpU->GetXaxis()->SetTitleFont(42);
  hFinalExpU->GetXaxis()->SetLabelFont(42);
  hFinalExpU->GetXaxis()->SetTitleSize(0.05);
  hFinalExpU->GetXaxis()->SetLabelSize(0.05);
  hFinalExpU->GetYaxis()->SetTitleFont(42);
  hFinalExpU->GetYaxis()->SetLabelFont(42);
  hFinalExpU->GetYaxis()->SetTitleSize(0.05);
  hFinalExpU->GetYaxis()->SetLabelSize(0.05);
  hFinalExpU->GetYaxis()->SetNdivisions(505);
  hFinalExpU->GetXaxis()->SetTitle("jet p_{T} (GeV)");
  hFinalExpU->GetYaxis()->SetTitle("(Data - Theory)/Theory");
  hFinalExpU->SetTitle("");
  for(i=0;i<hData->GetNbinsX();i++)
    {
      double a = hTotExpU->GetBinContent(i+1);
      double r = hData->GetBinContent(i+1);
      double e = r+a*(1+r);
      hFinalExpU->SetBinContent(i+1,e);
      a = hTotExpD->GetBinContent(i+1);
      r = hData->GetBinContent(i+1);
      e = r-fabs(a)*(r+1);
      hFinalExpD->SetBinContent(i+1,e);
    }
  hFinalExpU->SetFillColor(5);
  hFinalExpU->SetLineColor(5);
  hFinalExpU->SetLineWidth(1);
  hFinalExpD->SetFillColor(5);
  hFinalExpD->SetLineColor(5);
  hFinalExpD->SetLineWidth(1);
  hFinalTheoryU->SetFillColor(2);
  hFinalTheoryU->SetFillStyle(3013);
  hFinalTheoryU->SetLineColor(2);
  hFinalTheoryU->SetLineWidth(1);
  hFinalTheoryD->SetFillColor(2);
  hFinalTheoryD->SetFillStyle(3013);
  hFinalTheoryD->SetLineColor(2);
  hFinalTheoryD->SetLineWidth(1);
  hData->SetMarkerStyle(20);
  hData->SetLineWidth(1);
  hFinalExpU->Draw("AH HIST");
  hFinalExpD->Draw("same HIST");
  hFinalTheoryU->Draw("same HIST");
  hFinalTheoryD->Draw("same HIST");
  hFinalExpU->Draw("AXIS same HIST");
  TF1 *unit = new TF1("unit","0*x",0,2000);
  unit->SetLineColor(1);
  unit->SetLineStyle(2);
  unit->Draw("same");
  hData->Draw("P same");
  //TGraphAsymmErrors *gData = SetProperError(hData,10,dY);
  //gData->SetMarkerColor(4);
  //gData->Draw("P same Z");

  TPaveText *pave3 = new TPaveText(0.2,0.65,0.5,0.9,"NDC");
  pave3->SetTextFont(42);
  pave3->AddText("CMS preliminary");
  pave3->AddText(algoTitle);
  pave3->AddText(yTitle);
  pave3->AddText("#sqrt{s} = 10 TeV");
  pave3->SetFillColor(0);
  pave3->SetBorderSize(0);
  pave3->SetTextSize(0.04);
  pave3->Draw();

  TLegend *leg = new TLegend(0.51,0.65,0.94,0.9);
  leg->SetTextFont(42);
  leg->AddEntry(hData,"Data ( 10 pb^{-1} )","LP");
  leg->AddEntry(hFinalExpU,"Experimental Uncertainty","F");
  leg->AddEntry(hFinalTheoryU,"Theory Uncertainty","F");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  sprintf(canname,"FinalComparison_Y%d.eps",Ybin);
  cdata->Print(ALGO+"_"+canname);
}

void ClearHisto(TH1F* h, double Xmin, double Xmax)
{
  int i;
  double x;
  for(i=0; i<h->GetNbinsX();i++)
    {
      x = h->GetBinLowEdge(i+1);
      if (x<Xmin || x>=Xmax)
        {
          h->SetBinContent(i+1,0);
          h->SetBinError(i+1,0); 
        }
    } 
}

void AddQuadratic(TH1D *h1, TH1D *h2)
{
  int i;
  double y,y1,y2;
  for(i=0;i<h1->GetNbinsX();i++)
    {
      y1 = h1->GetBinContent(i+1);
      y2 = h2->GetBinContent(i+1);
      y = sqrt(y1*y1+y2*y2);
      h1->SetBinContent(i+1,y);
    }
}

TGraphAsymmErrors* SetProperError(TH1D *h, double Lumi, double dy)
{
  int i;
  double Nexp,dx,dy,nl,nh,a;
  double vx[200],vy[200],vexh[200],vexl[200],veyh[200],veyl[200];
  TGraphAsymmErrors *g;
  a = a = 0.3173/2;
  for(int i=0;i<h->GetNbinsX();i++)
    {
      vy[i] = h->GetBinContent(i+1);
      dx = h->GetBinWidth(i+1);
      Nexp = vy[i]*Lumi*dx*dy;
      vx[i] = h->GetBinCenter(i+1);
      vexh[i] = dx/2; 
      vexl[i] = dx/2;
      if (Nexp>20 || Nexp==0)
        {
          veyh[i] = h->GetBinError(i+1);
          veyl[i] = h->GetBinError(i+1);
        }
      else 
        {
          nl = Nexp-0.5*TMath::ChisquareQuantile(a,2*Nexp);
          nh = 0.5*TMath::ChisquareQuantile(1-a,2*(Nexp+1))-Nexp;
          veyl[i] = nl/(Lumi*dx*dy);
          veyh[i] = nh/(Lumi*dx*dy);
          cout<<Nexp<<" "<<veyl[i]<<" "<<veyh[i]<<endl;
        }
    }
  g = new TGraphAsymmErrors(i,vx,vy,vexl,vexh,veyl,veyh);
  return g; 
}


