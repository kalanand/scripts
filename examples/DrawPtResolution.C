
void DrawPtResolution(const char algorithm[100],const char njets[100])
{
  
  const double Pt[18]={5,15,25,35,45,60,75,90,120,150,200,300,400,550,750,1000,1500,5000};
  const int NETA = 6;
  const double eta_boundaries[NETA+1]={0.0,0.55,1.1,1.7,2.5,3.2,5.0};
  TFile *f;
  TH1F *hGenResolution[18],*hGenResp[18],*h,*hRef;
  TCanvas *c[NETA],*c1[NETA];
  int i,j,etabin,counter;
  char name[100],cname[100],filename[100];
  double r,e,x,ex;
  double vx[20],vy[20],vex[20],vey[20];
  TGraphErrors *gRes[NETA];
  TF1 *fitRes[NETA];

  sprintf(filename,"InclusiveNoteResults_%s_%sJets.root",algorithm,njets);
  f = new TFile(filename,"r");
  for(etabin=0;etabin<NETA;etabin++)
    { 
      sprintf(name,"ResolutionEtaBin%d",etabin);
      hGenResolution[etabin] = new TH1F(name,name,17,Pt);
      
      sprintf(name,"ResponseEtaBin%d",etabin);
      hGenResp[etabin] = new TH1F(name,name,17,Pt);  
    }

  for(i=0;i<17;i++)
   {
     sprintf(name,"Resolution_vs_Eta_RefPt%d",i);
     h = (TH1F*)f->Get(name);
     for(etabin=0;etabin<NETA;etabin++)
       {
         r = h->GetBinContent(etabin+1);
         e = h->GetBinError(etabin+1);
         hGenResolution[etabin]->SetBinContent(i+1,r);
         hGenResolution[etabin]->SetBinError(i+1,e);

       }
     
     sprintf(name,"Response_vs_Eta_RefPt%d",i);
     h = (TH1F*)f->Get(name);
     for(etabin=0;etabin<NETA;etabin++)
       {
         r = h->GetBinContent(etabin+1);
         e = h->GetBinError(etabin+1);
         hGenResp[etabin]->SetBinContent(i+1,r);
         hGenResp[etabin]->SetBinError(i+1,e);
       }
   }

for(etabin=0;etabin<NETA;etabin++)
{
  counter = 0;
  sprintf(name,"MeanRefPt_Eta%d",etabin); 
  hRef = (TH1F*)f->Get(name);  
  for(i=0;i<17;i++)
   {   
     sprintf(name,"Resolution_vs_Eta_RefPt%d",i);
     h = (TH1F*)f->Get(name);
     r = h->GetBinContent(etabin+1);
     e = h->GetBinError(etabin+1);
     x = hRef->GetBinContent(i+1);
     ex = hRef->GetBinError(i+1);
     if (r>0)
       {
         vy[counter] = r;
         vey[counter] = e;
         vx[counter] = x;
         vex[counter] = ex;
         counter++;
       } 
     gRes[etabin] = new TGraphErrors(counter,vx,vy,vex,vey);
     sprintf(name,"ResolutionFit_Eta%d",etabin); 
     fitRes[etabin] = new TF1(name,"sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",0,2000);
     fitRes[etabin]->SetParameters(0.04,1,1);
     fitRes[etabin]->SetParNames("C","S","N");
     fitRes[etabin]->SetLineColor(2);
     fitRes[etabin]->SetLineWidth(2);
     gRes[etabin]->Fit(fitRes[etabin],"RQ");
   } 
}
/*
for(etabin=0;etabin<NETA;etabin++)
        {
      sprintf(cname,"%s_PtClosure_EtaBin%d",algorithm,etabin);
      c[etabin] = new TCanvas(cname,cname,900,600);
      c[etabin]->cd();
      gPad->SetLogx();
      TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
      hGenResp[etabin]->SetTitle("Corrected jet response vs p_{T}^{gen}");
      hGenResp[etabin]->GetXaxis()->SetTitle("<p_{T}^{gen}> (GeV)");
      hGenResp[etabin]->GetYaxis()->SetTitle("<#frac{p_{T}^{cor}}{p_{T}^{gen}}>");
      hGenResp[etabin]->GetYaxis()->SetNdivisions(505);
      hGenResp[etabin]->SetMaximum(1.1);
      hGenResp[etabin]->SetMinimum(0.9);
      hGenResp[etabin]->SetMarkerStyle(21);
      hGenResp[etabin]->SetMarkerColor(1);
      hGenResp[etabin]->SetLineColor(1);
      hGenResp[etabin]->Draw();
      //Ratio->Draw("same");
      sprintf(name,"%1.1f<|#eta|<%1.1f",eta_boundaries[etabin],eta_boundaries[etabin+1]);
      leg->AddEntry(hGenResp[etabin],name,"LP");
      leg->SetFillColor(0);
      leg->SetHeader(algorithm);
      leg->Draw();   
}
*/
TPaveText *pave1 = new TPaveText(0.2,0.7,0.5,0.9,"NDC");
      pave1->AddText("Monte Carlo Truth");
      if (strcmp(algorithm,"SC7")==0)
        pave1->AddText("Seedless Cone R = 0.7");
      else
        pave1->AddText("k_{T} D = 0.6");
      pave1->SetFillColor(0);
      pave1->SetLineColor(0);
      pave1->SetBorderSize(0);
      pave1->SetTextFont(42);

TCanvas *closure = new TCanvas("Closure","Closure");
//gPad->SetLogx();
TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
hGenResp[0]->SetTitle("");
hGenResp[0]->GetXaxis()->SetTitle("Particle jet p_{T} (GeV)");
hGenResp[0]->GetXaxis()->SetRangeUser(100,1450);
hGenResp[0]->GetYaxis()->SetTitle("Corrected jet response");
hGenResp[0]->SetMaximum(1.04);
hGenResp[0]->SetMinimum(0.96);
hGenResp[0]->SetMarkerStyle(20);
hGenResp[0]->SetMarkerColor(1);
hGenResp[0]->SetLineColor(1);
hGenResp[2]->SetMarkerStyle(3);
hGenResp[2]->SetMarkerColor(1);
hGenResp[2]->SetLineColor(1);
hGenResp[3]->SetMarkerStyle(24);
hGenResp[3]->SetMarkerColor(1);
hGenResp[3]->SetLineColor(1);
hGenResp[0]->Draw();
hGenResp[2]->Draw("same");
hGenResp[3]->Draw("same");
sprintf(name,"%1.2f<|y|<%1.2f",eta_boundaries[0],eta_boundaries[1]);
leg->AddEntry(hGenResp[0],name,"LP");
sprintf(name,"%1.2f<|y|<%1.2f",eta_boundaries[2],eta_boundaries[3]);
leg->AddEntry(hGenResp[2],name,"LP");
sprintf(name,"%1.2f<|y|<%1.2f",eta_boundaries[3],eta_boundaries[4]);
leg->AddEntry(hGenResp[3],name,"LP");
leg->SetFillColor(0);
leg->SetBorderSize(0); 
leg->SetTextFont(42);
leg->Draw();   
pave1->Draw();
TF1 *fit1 = new TF1("fit1","0*x+1",0,2000);
fit1->SetLineStyle(2);
fit1->SetLineColor(1);
fit1->Draw("same");
/*
  for(etabin=0;etabin<NETA;etabin++)
    {
      sprintf(cname,"%s_Resolution_EtaBin%d",algorithm,etabin);
      c1[etabin] = new TCanvas(cname,cname,900,600);
      c1[etabin]->cd();
      gPad->SetLogx();
      sprintf(name,"%1.1f<#eta<%1.1f",eta_boundaries[etabin],eta_boundaries[etabin+1]);
      hGenResolution[etabin]->SetTitle(name);
      hGenResolution[etabin]->GetXaxis()->SetTitle("<p_{T}^{gen}> (GeV)");
      hGenResolution[etabin]->GetYaxis()->SetTitleOffset(1.1);   
      hGenResolution[etabin]->GetYaxis()->SetTitle("#sigma(p_{T}^{cor}/p_{T}^{gen})");
      hGenResolution[etabin]->GetYaxis()->SetNdivisions(505); 
      hGenResolution[etabin]->SetMarkerStyle(20);
      hGenResolution[etabin]->SetMarkerColor(1);
      hGenResolution[etabin]->SetLineColor(1);
      hGenResolution[etabin]->Draw("same");
    }

sprintf(filename,"ResolutionResults_%s_%sJets.root",algorithm,njets);
TFile *outf = new TFile(filename,"RECREATE");
outf->cd();
for(etabin=0;etabin<NETA;etabin++)
    {
      sprintf(cname,"%s_Resolution_EtaBin%d",algorithm,etabin);
      c1[etabin] = new TCanvas(cname,cname,900,600);
      c1[etabin]->cd();
      gPad->SetLogx();
      sprintf(name,"%1.1f<#eta<%1.1f",eta_boundaries[etabin],eta_boundaries[etabin+1]);
      gRes[etabin]->SetTitle(name);
      gRes[etabin]->GetXaxis()->SetTitle("<p_{T}^{gen}> (GeV)");
      gRes[etabin]->GetYaxis()->SetTitleOffset(1.1);   
      gRes[etabin]->GetYaxis()->SetTitle("#sigma(p_{T}^{cor}/p_{T}^{gen})");
      gRes[etabin]->GetYaxis()->SetNdivisions(505); 
      gRes[etabin]->SetMarkerStyle(20);
      gRes[etabin]->SetMarkerColor(1);
      gRes[etabin]->SetLineColor(1);
      gRes[etabin]->Draw("AP");
      sprintf(name,"ResolutionGraph_Eta%d",etabin);
      gRes[etabin]->Write(name); 
      fitRes[etabin]->Write();  
    }
*/
}
