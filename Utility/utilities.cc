
const double M_PI = TMath::Pi();



Float_t CosThetaStar(TLorentzVector *daughter, TLorentzVector *parent)
{
  // angle between mother and daughter in mother's rest frame
  
  TVector3 boostvect = -(*parent).BoostVector();

  //boost the daughter to mother's rest frame
  (*daughter).Boost(boostvect);

  // get the unit vectors for the mother and daughter
  TVector3 dunit = (*daughter).Vect().Unit();
  TVector3 punit = (*parent).Vect().Unit();

  // take their dot product
  return dunit.Dot(punit);
}






// calculate m_q* using the (pT)_jet = (pT)_Z
TLorentzVector* qstarMassFromPtBalance(TLorentzVector *Z, 
				       double Jeteta){

  TLorentzVector jet( -Z->Px(), -Z->Py(), Z->Pt()*sinh(Jeteta), 
		      Z->Pt()*cosh(Jeteta) );

  TLorentzVector* qstar = new TLorentzVector(*Z);
  *qstar += jet;
  return qstar;
}




// Create a user-defined histogram
TH1D* CreateHistogram(TString name, int bin, double min, 
		      double max, int color=1, 
		      TString xtit="", TString ytit="", bool sumw2=true) {

  TH1D* h = new TH1D(name,"", bin, min, max);
  if(sumw2) h->Sumw2();
  h->SetLineWidth(2);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  SetTitles( *h, xtit, ytit);

  return h;
}







// Set histogram axis titles and margins
void SetTitles(TH1& h, TString xtit, TString ytit) {

  TAxis* xax = h.GetXaxis();
  TAxis* yax = h.GetYaxis();
  xax->SetTitle(xtit);
  xax->SetTitleSize(0.04);
  xax->SetTitleOffset(1.4);
  xax->SetNdivisions(505);
  yax->SetTitle(ytit);
  yax->SetTitleOffset(1.3);
}






// Plot on canvas: one histogram
TCanvas* PlotOnCanvas(TString plotname, TH1& h, char* drawOption="e") {
  TCanvas* can = new TCanvas( plotname, plotname, 500, 500 );
  gStyle->SetOptStat(0);  
  h.Draw(drawOption);
  can->SaveAs( plotname+TString(".eps") ); 
  can->SaveAs( plotname+TString(".gif") );  
  can->SaveAs( plotname+TString(".root") ); 
  return can; 
}





// Plot on canvas: two histograms
TCanvas* PlotTwoOnCanvas(TString plotname, TH1& h1, TString label1,
			 TH1& h2, TString label2, float x0=0.63, 
			 float y0=0.8, float x1=0.9, float y1=0.9, 
			 char* opt1 = "HIST e", char* opt2 = "Hist esame") {
  
  double i1 = h1.Integral();
  double i2 = h2.Integral();
  
  TCanvas* can = new TCanvas( plotname, plotname, 500, 500 );    
  gStyle->SetOptStat(0); 
  if( !(i1==0.0 || i2==0.0) ) {
    h1.Draw("Hist e");  
    h2.Draw("Hist esame"); 
  }
  else if(i2==0.0) h1.Draw("e");
  else h2.Draw("e");
  TLegend* leg = new TLegend(x0,y0,x1,y1);
  leg->AddEntry( &h1, label1,"L"); 
  leg->AddEntry( &h2, label2,"L"); 
  leg->SetMargin(0.15);
  leg->SetFillColor(0);
  if( !(i1==0.0 || i2==0.0) ) leg->Draw();
  // gPad->SetLogy();
 //  gPad->SetGridx();  
  can->SaveAs( plotname+TString(".eps") ); 
  can->SaveAs( plotname+TString(".gif") );  
  can->SaveAs( plotname+TString(".root") );
  return can; 
}





// Plot on canvas: three histograms
TCanvas* PlotThreeOnCanvas(TString plotname, TH1& h1, TString label1,
			   TH1& h2, TString label2, TH1& h3, TString label3,
			   float x0=0.63, 
			   float y0=0.8, float x1=0.9, float y1=0.9) {
  
  TCanvas* can = new TCanvas( plotname, plotname, 500, 500 );
  h1.Draw("Hist e");  
  h2.Draw("Hist esame"); 
  h3.Draw("Hist esame");
  TLegend* leg = new TLegend(x0,y0,x1,y1);
  leg->AddEntry( &h1, label1,"L"); 
  leg->AddEntry( &h2, label2,"L"); 
  leg->AddEntry( &h3, label3,"L"); 
  leg->SetMargin(0.15);
  leg->SetFillColor(0);
  leg->Draw();
  // gPad->SetLogy();
  gPad->SetGridx();  
  can->SaveAs( plotname+TString(".eps") ); 
  can->SaveAs( plotname+TString(".gif") );  
  can->SaveAs( plotname+TString(".root") );
  return can; 
}




// Plot on canvas: four histograms
TCanvas* PlotFourOnCanvas(TString plotname, TH1& h1, TString label1,
			  TH1& h2, TString label2, TH1& h3, TString label3,
			  TH1& h4, TString label4, float x0=0.63, 
			  float y0=0.8, float x1=0.9, float y1=0.9) {
  
  TCanvas* can = new TCanvas( plotname, plotname, 500, 500 );
  h1.Draw("Hist e");  
  h2.Draw("Hist esame"); 
  h3.Draw("Hist esame");
  h4.Draw("Hist esame");
  TLegend* leg = new TLegend(x0,y0,x1,y1);
  leg->AddEntry( &h1, label1,"L"); 
  leg->AddEntry( &h2, label2,"L"); 
  leg->AddEntry( &h3, label3,"L"); 
  leg->AddEntry( &h4, label4,"L"); 
  leg->SetMargin(0.15);
  leg->SetFillColor(0);
  leg->Draw();
  // gPad->SetLogy();
  //   gPad->SetGridx(); 
  
  can->SaveAs( plotname+TString(".eps") ); 
  can->SaveAs( plotname+TString(".gif") );  
  can->SaveAs( plotname+TString(".root") );
  return can; 
}





// Plot on canvas: five histograms
TCanvas* PlotFiveOnCanvas(TString plotname, TH1& h1, TString label1,
			  TH1& h2, TString label2, TH1& h3, TString label3,
			  TH1& h4, TString label4, TH1& h5, TString label5, 
			  float x0=0.63, 
			  float y0=0.8, float x1=0.9, float y1=0.9) {
  
  TCanvas* can = new TCanvas( plotname, plotname, 500, 500 );
  h1.Draw("Hist e");  
  h2.Draw("Hist esame"); 
  h3.Draw("Hist esame");
  h4.Draw("Hist esame");
  h5.Draw("Hist esame");
  TLegend* leg = new TLegend(x0,y0,x1,y1);
  leg->AddEntry( &h1, label1,"L"); 
  leg->AddEntry( &h2, label2,"L"); 
  leg->AddEntry( &h3, label3,"L"); 
  leg->AddEntry( &h4, label4,"L"); 
  leg->AddEntry( &h5, label5,"L"); 
  leg->SetMargin(0.15);
  leg->SetFillColor(0);
  leg->Draw();
  // gPad->SetLogy();
  gPad->SetGridx();  
  can->SaveAs( plotname+TString(".eps") ); 
  can->SaveAs( plotname+TString(".gif") );  
  can->SaveAs( plotname+TString(".root") );
  return can; 
}




// find the leading and second jet indices

void FindLeadIndex(float pT[], float eta[], float phi[],
		   float e1eta, float e1phi, 
		   float e2eta, float e2phi, int &lead, int &sec) {
  float max = 0.0;
  lead = -1;
  for (int i=0; i<10; i++) {
    
    double r1 = radius( eta[i],phi[i], e1eta, e1phi );
    double r2 = radius( eta[i],phi[i], e2eta, e2phi );
    if(!(r1>0.2 && r2>0.2) ) continue;
    if(pT[i]>max) { max = pT[i]; lead = i; }
  }

  max = 0.0;
  sec = -1;
  for (int i=0; i<10; i++) {
    if(i==lead) continue;
    double r1 = radius( eta[i],phi[i], e1eta, e1phi );
    double r2 = radius( eta[i],phi[i], e2eta, e2phi );
    if(!(r1>0.2 && r2>0.2) ) continue;
    if(pT[i]>max) { max = pT[i]; sec = i; }
  }
}





// find the leading and second jet indices

void FindLeadIndex(float pT[], int &lead, int &sec) {

  float max = 0.0;
  lead = -1;
  for (int i=0; i<10; i++) {
    if(pT[i]>max) { max = pT[i]; lead = i; }
  }

  max = 0.0;
  sec = -1;
  for (int i=0; i<10; i++) {
    if(i==lead) continue;
    if(pT[i]>max) { max = pT[i]; sec = i; }
  }
}





