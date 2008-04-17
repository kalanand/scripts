// study the NRPW

void NRPW() {
  BdkPdfDDalitz w("w", "w", *m12, *m13, BdkDalitzBase::D0, 0, 
		  BdkDDalitzAmp::NRPW_PM + 
		  BdkDDalitzAmp::NRPW_M0 + 
		  BdkDDalitzAmp::NRPW_0P);

  TCanvas * can = new TCanvas("can", "can", 1000, 500);
  can->Divide(2,1);

  /*
  w.pdfType()->dalitzAmp()->calNorm();

  w.parameters().readFromFile("NRPW.par");
  w.parameters().Print("V");

  w.pdfType()->dalitzAmp()->fitFractions().Print("V");

  data = w.generate(10000); 

  TH2 * hist = data->createHistogram(*m12, *m13);
  can->cd(1);
  hist->Draw();
*/


  BdkPdfDDalitz r("r", "r", *m12, *m13, BdkDalitzBase::D0, 0, 
		  BdkDDalitzAmp::NRPW_PM + 
		  BdkDDalitzAmp::NRPW_M0 + 
		  BdkDDalitzAmp::NRPW_0P + 
		  BdkDDalitzAmp::RHO0 + 
		  BdkDDalitzAmp::RHOP + 
		  BdkDDalitzAmp::RHOM);

  r.pdfType()->dalitzAmp()->calNorm();

  r.parameters().readFromFile("NRPW.par");
  r.parameters().Print("V");

  r.pdfType()->dalitzAmp()->fitFractions().Print("V");

  RooDataSet * d2 = r.generate(10000); 
  TH2 * hist2 = d2->createHistogram(*m12, *m13);
  can->cd(2);
  hist2->Draw();
}
  
/*
NRPW.par contains the following:
w.pdf.dalitzAmp.NRPW_0P_amp =  15.626 +/- 1.7199 
w.pdf.dalitzAmp.NRPW_0P_phase = -51.3850 +/- 7.4643 
w.pdf.dalitzAmp.NRPW_M0_amp =  16.298 +/- 1.6898 
w.pdf.dalitzAmp.NRPW_M0_phase = -51.8422 +/- 7.1470 
w.pdf.dalitzAmp.NRPW_PM_amp =  15.784 +/- 1.6688 
w.pdf.dalitzAmp.NRPW_PM_phase = -50.6650 +/- 7.4565 
*/
