
void DrawPtResolutionComparison(const char algorithm[100])
{
  
  
  const int NETA = 6;
  const double eta_boundaries[NETA+1]={0.0,0.55,1.1,1.7,2.5,3.2,5.0};
    
  TCanvas *c[NETA],*c1[NETA];
  int i,j,etabin,counter;
  char name[100],cname[100],filename[100];
  
  TF1 *fitRes[NETA][2];

  TFile *inf[2];
  sprintf(filename,"ResolutionResults_%s_AllJets.root",algorithm);
  inf[0] = new TFile(filename);
  sprintf(filename,"ResolutionResults_%s_LeadingJets.root",algorithm);
  inf[1] = new TFile(filename);

  for(etabin=0;etabin<NETA;etabin++)
    { 
      sprintf(name,"ResolutionFit_Eta%d",etabin); 
      fitRes[etabin][0] = (TF1*)inf[0]->Get(name);
      fitRes[etabin][1] = (TF1*)inf[1]->Get(name);
    }

for(etabin=0;etabin<NETA;etabin++)
    {
      double ptmax = TMath::Min(2000.,5000./cosh(eta_boundaries[etabin]));
      TH1F *tmp = new TH1F("tmp","tmp",100,20,ptmax);
      sprintf(name,"%1.2f<|y|<%1.2f",eta_boundaries[etabin],eta_boundaries[etabin+1]); 
      TPaveText *pave = new TPaveText(0.6,0.7,0.9,0.9,"NDC");
      pave->AddText("Monte Carlo Truth");
      pave->AddText(name);
      if (strcmp(algorithm,"SC7")==0)
        pave->AddText("Seedless Cone R = 0.7");
      else
        pave->AddText("k_{T} D = 0.6");
      pave->SetFillColor(0);
      pave->SetLineColor(0);
      pave->SetBorderSize(0);
      pave->SetTextFont(42);
      TLegend *leg = new TLegend(0.6,0.6,0.9,0.7);

      sprintf(cname,"%s_Resolution_EtaBin%d",algorithm,etabin);
      c1[etabin] = new TCanvas(cname,cname);
      c1[etabin]->cd();
      gPad->SetLogx();
      tmp->SetTitle("");
      tmp->SetMaximum(0.4);
      tmp->SetMinimum(0.);  
      tmp->GetXaxis()->SetTitle("Particle jet p_{T} (GeV)");
      tmp->GetYaxis()->SetTitle("Relative Energy Resolution");
      fitRes[etabin][0]->SetLineColor(1);
      fitRes[etabin][1]->SetLineColor(2);
      tmp->Draw(); 
      fitRes[etabin][0]->Draw("same");
      fitRes[etabin][1]->Draw("same");
      leg->AddEntry(fitRes[etabin][0],"All jets","L");
      leg->AddEntry(fitRes[etabin][1],"2 leading jets","L");
      leg->SetFillColor(0);
      leg->SetLineColor(0);
      leg->SetBorderSize(0); 
      leg->SetTextFont(42);
      leg->Draw();
      pave->Draw();
    }

}
