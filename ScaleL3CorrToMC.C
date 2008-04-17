void ScaleL3CorrToMC () {

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(000); 
  gStyle->SetPalette(1);
 

  double x,y,e;

  double rMin = 5.0;
  double rMax = 5000.0;
  int rBin = (int)(rMax - rMin);

  TFile* inf = new TFile("L3Graphs_test_Icone5.root","r");
  TGraphErrors *g_Cor = (TGraphErrors*)inf->Get("Correction_vs_CaloPt");
  TMatrixD *COV_Cor = (TMatrixD*)inf->Get("CovMatrix_Correction");
  TF1 *CorFit = (TF1*)g_Cor->GetFunction("CorFit");
  CorFit->SetRange(rMin,rMax);  


 TF1* MCcurve = new TF1("MCcurve","[0]+[1]/(pow(log10(x),[2])+[3])", 
			rMin, rMax); 
 MCcurve->SetParameter(0, 0.996998);
 MCcurve->SetParameter(1, 4.39412);
 MCcurve->SetParameter(2, 2.96134);
 MCcurve->SetParameter(3, 1.69966);
 MCcurve->SetLineColor(4);
 MCcurve->SetLineStyle(1);


 TF1* MCAsData = new TF1("MCAsData","[0]+[1]/(pow(log10(x),[2])+[3])", 
			rMin, rMax); 
 MCAsData->SetParameter(0, CorFit->GetParameter(0));
 MCAsData->SetParameter(1, CorFit->GetParameter(1));
 MCAsData->SetParameter(2, CorFit->GetParameter(2));
 MCAsData->SetParameter(3, CorFit->GetParameter(3));


 TF1* MCOverCor = new TF1("MCOverCor","MCcurve/MCAsData",rMin, rMax); 
 MCOverCor->SetLineColor(6);
 MCOverCor->SetLineStyle(1);
 MCOverCor->SetLineWidth(2);



 /*
 TF1* closure = new TF1("closure","MCAsData*MCOverCor",rMin, rMax); 
 closure->SetLineColor(3);
 closure->SetLineStyle(2);
 closure->SetLineWidth(2);
*/



  hCorUncertainty = new TH1F("CorrectionUncertainty",
			     "CorrectionUncertainty",
			     rBin, rMin, rMax);
  for(int i=0;i<rBin;i++) {
      x = hCorUncertainty->GetBinCenter(i+1);
      y = CorFit->Eval(x);
      e = FitUncertainty(CorFit,COV_Cor,x);
      hCorUncertainty->SetBinContent(i+1,y);
      hCorUncertainty->SetBinError(i+1,e); 
    }




  //////////////////////////////////////////////////
  const Int_t nPoints = 16;
  Double_t new_gr_x[nPoints];
  Double_t new_gr_y[nPoints];
  Double_t new_gr_errx[nPoints];
  Double_t new_gr_erry[nPoints];
  Double_t grx[nPoints] = g_Cor->GetX();


  for(int i=0; i<nPoints; ++i) {
    new_gr_x[i] = g_Cor->GetX()[i];
    new_gr_y[i] = MCcurve->Eval(new_gr_x[i]);
    new_gr_errx[i] = g_Cor->GetErrorX(i); 
    new_gr_erry[i] = g_Cor->GetErrorY(i)*(MCOverCor->Eval(new_gr_x[i])); 
  }

  TGraphErrors* new_gr = 
    new TGraphErrors(nPoints, new_gr_x, new_gr_y, 
		     new_gr_errx, new_gr_erry);



  newhCorUncertainty = new TH1F("newCorrectionUncertainty",
			     "", rBin, rMin, rMax);
  for(int i=0;i<rBin;i++) {
      x = newhCorUncertainty->GetBinCenter(i+1);
      y = CorFit->Eval(x)*(MCOverCor->Eval(x));
      e = FitUncertainty(CorFit,COV_Cor,x);
      newhCorUncertainty->SetBinContent(i+1,y);
      newhCorUncertainty->SetBinError(i+1,e); 
    }

  //////////////////////////////////////////////////


  TCanvas *c_Correction = new TCanvas("Correction","Correction",900,600);
  c_Correction->cd(); 
  gPad->SetLogx();
  hCorUncertainty->SetTitle(""); 
  g_Cor->SetMarkerStyle(20);
  g_Cor->SetMarkerColor(1);
  g_Cor->SetLineColor(1);
  hCorUncertainty->SetMaximum(3.2);
  hCorUncertainty->SetMinimum(1);
  hCorUncertainty->GetYaxis()->SetNdivisions(505);
  hCorUncertainty->GetXaxis()->SetTitle("p_{T} (GeV)");
  hCorUncertainty->GetYaxis()->SetTitle("L3Correction factor");
  CorFit->SetLineColor(2);
  CorFit->SetLineWidth(2);
  hCorUncertainty->SetLineColor(5);
  hCorUncertainty->SetFillColor(5);
  hCorUncertainty->SetMarkerColor(5);
  newhCorUncertainty->SetLineColor(5);
  newhCorUncertainty->SetFillColor(5);
  newhCorUncertainty->SetMarkerColor(5);


  const int  plotCase=2;

  if(plotCase==1) {
  // hCorUncertainty->Draw("E3");
  newhCorUncertainty->Draw("E3");

  //g_Cor->Draw("Psame"); 
  new_gr->Draw("Psame"); 
  MCcurve->Draw("Lsame");
  //  MCOverCor->Draw("Lsame");
  //  closure->Draw("Lsame");
  TLegend *leg = new TLegend(0.5,0.65,0.89,0.89);
  leg->AddEntry(g_Cor,"L3 measurement from Z#rightarrow e^{+}e^{-}+Jet","LP");
  //  leg->AddEntry(CorFit,"fit","L");
  leg->AddEntry(MCcurve,"fit","L");
  // leg->AddEntry(hCorUncertainty,"fit uncertainty","F");
  leg->AddEntry(newhCorUncertainty,"fit uncertainty","F");
  //leg->AddEntry(MCcurve,"L3Correction from dijet MC-truth","L");
  // leg->AddEntry(MCOverCor,"Flavor-map function","L");
  // leg->AddEntry(closure,"Flavor-map function #times fit","L");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
  }

  if(plotCase==2) {
    newhCorUncertainty->Draw("AXIS");
    MCOverCor->Draw("Lsame"); 
    TLegend *leg = new TLegend(0.5,0.65,0.89,0.89);
    leg->AddEntry(MCOverCor,"Flavor-map function","L");
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->Draw();
  }

  if(plotCase==3) {
    new_gr->Draw("AP"); 
    MCcurve->Draw("Lsame");
    TLegend *leg = new TLegend(0.5,0.65,0.89,0.89);
    leg->AddEntry(g_Cor,"L3 measurement from Z#rightarrow e^{+}e^{-}+Jet","LP");
    leg->AddEntry(MCcurve,"L3Correction from dijet MC-truth","L");
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->Draw();
  }


  if(plotCase==4) {
    newhCorUncertainty->Draw("E3");
    new_gr->Draw("Psame"); 
    new_gr->Fit(MCcurve,"RQ");
    MCcurve->SetLineColor(2);
    MCcurve->Draw("Lsame");
    TLegend *leg = new TLegend(0.5,0.65,0.89,0.89);
    leg->AddEntry(g_Cor,"L3 measurement from Z#rightarrow e^{+}e^{-}+Jet","LP");
    leg->AddEntry(MCcurve,"fit","L");
    leg->AddEntry(newhCorUncertainty,"fit uncertainty","F");
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->Draw();
  }


}




double FitUncertainty( TF1* f, TMatrixD* COV, double x)
{
  int i,j;
  double df,sum,y,z,x;
  double PartialDerivative[10],Parameter[10];
  int npar = 4;
  int N = f->GetNumberFreeParameters();
  int dim = COV->GetNrows();

  if (dim != npar || N != npar)
    {
      cout<<"ERROR: wrong number of parameters !!!!"<<endl;
      return(-1);
    }  
  for(i=0;i<npar;i++)
    Parameter[i] = f->GetParameter(i);
  z = pow(log10(x),Parameter[2]);  
  PartialDerivative[0] = 1.;
  PartialDerivative[1] = 1./(z+Parameter[3]);
  PartialDerivative[3] = -Parameter[1]/pow(z+Parameter[3],2);
  PartialDerivative[2] = PartialDerivative[3]*log(log10(x))*z;
  sum = 0.;
  for(i=0;i<npar;i++)
    for(j=0;j<npar;j++)
      {
        y = PartialDerivative[i]*PartialDerivative[j]*COV(i,j);
        sum+=y;
      }
  df = sqrt(sum);
  return df;
}
