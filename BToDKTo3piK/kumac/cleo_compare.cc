// This file compares our D Dalitz PDF to the CLEO results. See end of file
// for the pars file it needs.
{
RooRealVar m12("m12", "m12", 0, 3);
RooRealVar m13("m13", "m13", 0, 3);

RooArgList list(m12,m13);

RooFormulaVar m23("m23", "m23", "1.8645*1.8645+0.1349766*0.1349766+0.13957*0.13957+0.13957*0.13957-m12-m13", list);


 BdkPdf2DpolyDalitz eff("eff", "eff", m12, m13);
 BdkDalitzBase::setEfficiencyFunc(&eff);
 eff.parameters().Print("V");


 BdkPdfDDalitz amp("amp", "amp", m12, m13, BdkDalitzBase::D0);

//amp.parameters().readFromFile("pars-cleo-flipped");
//amp.parameters().readFromFile("pars-cleo");
amp.parameters().readFromFile("pars-equal");

//amp.pdfType()->dalitzAmp()->calNorm(1000000);
//amp.parameters().Print("V");


RooDataSet * data = amp.generate(20000);
data->addColumn(m23);

ofstream file("data.dat");
for (int i = 0; i < data->numEntries(); ++i){
  file << ((RooAbsReal*)(data->get(i)->find("m12")))->getVal() << " "
       << ((RooAbsReal*)(data->get(i)->find("m13")))->getVal() << " "
       << ((RooAbsReal*)(data->get(i)->find("m23")))->getVal() << endl;
}

// Below, make 2D plots of the real, imaginary, phase, and amp^2 of the PDF:
/*
TH2F * hReal = new TH2F("hReal", "hReal", 
		       30, m12.getMin(), m12.getMax(), 
		       30, m13.getMin(), m13.getMax());
		      
TH2F * hImag = new TH2F("hImag", "hImag", 
		      30, m12.getMin(), m12.getMax(), 
		      30, m13.getMin(), m13.getMax());
		      
TH2F * hPhase = new TH2F("hPhase", "hPhase", 
		      30, m12.getMin(), m12.getMax(), 
		      30, m13.getMin(), m13.getMax());
		      
TH2F * hAmp2 = new TH2F("hAmp2", "hAmp2", 
		      30, m12.getMin(), m12.getMax(), 
		      30, m13.getMin(), m13.getMax());		      

const double step = 0.1;
for (double x = m12.getMin(); x < m12.getMax(); x += step) {
  for (double y = m13.getMin(); y < m13.getMax(); y += step) {
    if (amp.pdfType()->inDalitz(x,y)) {
      RooComplex comp = amp.pdfType()->dalitzAmp()->getamp(x, y);
      double real = comp.re();
      double imag = comp.im();
      double phase = atan2(imag, real);
      double amp2 = comp.abs2();
      
      hReal->Fill(x, y, real);
      hImag->Fill(x, y, imag);
      hPhase->Fill(x, y, phase);
      hAmp2->Fill(x, y, amp2);
    }
  }
}
	 
TCanvas * canvas = new TCanvas("can0", "can0", 1200, 1200);
canvas->Divide(2,2);
canvas->cd(1);
hReal->Draw("contz");
canvas->cd(2);
hImag->Draw("contz");
canvas->cd(3);
hPhase->Draw("contz");
canvas->cd(4);
hAmp2->Draw("contz");
*/
}

/*
Copy the stuff below to pars-cleo-flipped:

amp.pdf.dalitzAmp.Nonres_amp     =  1.0300 L(0 - 100)
amp.pdf.dalitzAmp.Nonres_phase   =  77.000 L(-10000 - 10000)
amp.pdf.dalitzAmp.Rho+_amp       =  1.0000 C
amp.pdf.dalitzAmp.Rho+_phase     =  0.0000 C
amp.pdf.dalitzAmp.Rho-_amp       =  0.65000 L(0 - 100)
amp.pdf.dalitzAmp.Rho-_phase     =  -4.00 L(-10000 - 10000)
amp.pdf.dalitzAmp.Rho0_amp       =  0.56000 L(0 - 100)
amp.pdf.dalitzAmp.Rho0_phase     =  10.000 L(-10000 - 10000)

To make pars-cleo, change the Rho-_phase from -4 to 176.

To make pars-equal, set the amplitudes of the 3 resonances to 1,
that of the non-resonant component to 0, and all the phases to 0.
*/
