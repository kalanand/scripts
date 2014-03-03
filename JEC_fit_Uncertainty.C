void JEC_fit_Uncertainty(int N)
{
  TF1 *func[1000];
  TFile *inf = new TFile("L3Graphs_test_Icone5.root");
  TGraphErrors *g = (TGraphErrors*)inf->Get("Correction_vs_CaloPt");
  TGraphErrors *vg[1000];
  int i,k;
  double x[20],y[20],ex[20],ey[20];
  double vx[20],vy[20],vex[20],vey[20];
  for(i=0;i<g->GetN();i++)
    {
      g->GetPoint(i,x[i],y[i]);
      ex[i]=g->GetErrorX(i);
      ey[i]=g->GetErrorY(i); 
    }  
  TRandom *rnd = new TRandom();
  rnd->SetSeed(0);
  for(k=0;k<N;k++)
    {
      for(i=0;i<g->GetN();i++)
        {	
          vx[i] = rnd->Gaus(x[i],ex[i]);
          //vx[i] = x[i];
          vy[i] = rnd->Gaus(y[i],ey[i]);
          vex[i] = ex[i];
          vey[i] = ey[i];
        }
      vg[k] = new TGraphErrors(g->GetN(),vx,vy,vex,vey);
      func[k] = new TF1("func","[0]+[1]/(pow(log10(x),[2])+[3])",1,2000);
      func[k]->SetParameters(1,3,6,5);
      vg[k]->Fit(func[k],"RQ");     	
    }
  
  TCanvas *c = new TCanvas("c","c");
  gPad->SetLogx();
  g->SetMarkerStyle(20);
  g->SetMaximum(3.5);
  g->Draw("AP");
  for(k=0;k<N;k++)
    {
      func[k]->SetLineColor(5);
      func[k]->SetLineWidth(1);
      cout<<func[k]->GetChisquare()<<endl;
      vg[k]->SetMarkerColor(2);
      vg[k]->SetLineColor(2);
      vg[k]->SetMarkerStyle(21);
      //if (func[k]->GetChisquare()<0.1)
        //vg[k]->Draw("sameP");
      func[k]->Draw("same");  
    }  	 
}  


