// fits the parameters of the BBBadD mES-mD 2D pdf

void getBBBadDMesMd(const char * parFile = "params/mes-d0-dep.par") {
  BdkPdf2CBACBP & pdf = pdfOnResDK.BBBadD0Prod().mESmDPdf();
  pdf.parameters().readFromFile(parFile);

  //  pdf.defaultLinks();
  pdf.linkCBAEndPoint();
  pdf.linkCBAM0();
  //pdf.linkCBAAlpha();
  pdf.linkCBAEnne();
  pdf.linkCBASigma();

  pdf.linkCBPAlpha();
  pdf.linkCBPEnne();
  //pdf.linkCBPM0();
  //pdf.linkCBPSigma();



  pdf.parameters().Print("V");

  readCut = cutBadD;
  data = read(chainBBComb);

  fitOption = "mrv";
  doFit = kFALSE;
  fit(pdf, *data);

  if (doPlot) {
    setupPlots();
    data->plotOn(mesFrame);
    pdf.getPdf()->plotOn(mesFrame);
    data->plotOn(d0massFrame);
    pdf.getPdf()->plotOn(d0massFrame);

    TCanvas * can = new TCanvas("can", "can", 1100, 500);
    can->Divide(2,1);
    can->cd(1);
    mesFrame->Draw();
    can->cd(2);
    d0massFrame->Draw();

    cout << "chi^2 of mes plot = " << mesFrame->chiSquare() << endl;
    cout << "chi^2 of mD plot  = " << d0massFrame->chiSquare() << endl;
  }
}
