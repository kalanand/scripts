// Function getPdfInfo provides information about the Pdf shapes of the 
// estimators for rho and theta and their errors as a function of the 
// true value of rho. It also tests whether these depend on the true value
// of theta (the dependence is very small). This provides information
// with which to calculate the confidence intervals.
//


const int nBinsPulls = 20;
const char * pullVsBinFileName = "pullVsBin.root";

void getPdfInfo(int nCharges = 1, // 2 for both B- and B+ 
		double x0 = 0.85, 
		double genSpread = 0.4) {

  // This is my link to the file /nfs/farm/babar/AWG17/BCK/Frank/toyMC/rndToy-2/fitMC.root:
  TFile * file = new TFile("fitMC-correct-for-pull-full-fit.root");
  RooDataSet * data = (RooDataSet*)file->Get("fitParData");
  TTree & tree = data->tree();

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);

  // The rho bin edges are determined by x0 +/- the spread
  // with which we generated the events:
  const double LOR = x0 - genSpread;
  const double HIR = x0 + genSpread;

  // plot correlations between fit and true variables:
  can = new TCanvas("corr", "corr", 800, 400 * nCharges); 
  can->Divide(2, nCharges);

  can->cd(1);
  tree.Draw("dalitzHolderN.sigGoodD0.rho:rhoNgen");

  can->cd(2);
  tree.Draw("dalitzHolderN.sigGoodD0.theta:thetaNgen");

  if (2 == nCharges) {
    can->cd(3);
    tree.Draw("dalitzHolderP.sigGoodD0.rho:rhoPgen");
    
    can->cd(4);
    tree.Draw("dalitzHolderP.sigGoodD0.theta:thetaPgen");
  }

  // Plot the correlation between rho and theta:
  can = new TCanvas("rhoethetacorr", "rhoe-theta corr", 400 * nCharges, 400); 
  can->Divide(nCharges, 1);
  can->cd(1);
  tree.Draw("(dalitzHolderN.sigGoodD0.rho-rhoNgen):(dalitzHolderN.sigGoodD0.theta-thetaNgen)", "abs(dalitzHolderN.sigGoodD0.theta-thetaNgen)/dalitzHolderN.sigGoodD0.thetaerr<7&&dalitzHolderN.sigGoodD0.thetaerr<90");

  if (2 == nCharges) {  
    can->cd(2);
    tree.Draw("(dalitzHolderP.sigGoodD0.rho-rhoPgen):(dalitzHolderP.sigGoodD0.theta-thetaPgen)", "abs(dalitzHolderP.sigGoodD0.theta-thetaPgen)/dalitzHolderP.sigGoodD0.thetaerr<7&&dalitzHolderP.sigGoodD0.thetaerr<90");
  }

  can->SaveAs("rho-theta-corr.eps");
    
  // profile error dependences on rho showing spread:
  can = new TCanvas("errProfRhoS", "Errors vs. true rho, bars=spreads", 800, 400 * nCharges); 
  can->Divide(2, nCharges);

  can->cd(1);
  TProfile * pRNErrR = new TProfile("pRNErrR", "RN err vs. RNgen", 50, LOR, HIR, "s");
  tree.Project("pRNErrR", "dalitzHolderN.sigGoodD0.rhoerr:rhoNgen", "dalitzHolderN.sigGoodD0.rhoerr<0.4");
  pRNErrR->Draw();

  can->cd(2);
  TProfile * pTNErrR = new TProfile("pTNErrR", "TN err vs. RNgen", 50, LOR, HIR, "s");
  tree.Project("pTNErrR", "dalitzHolderN.sigGoodD0.thetaerr:rhoNgen", "dalitzHolderN.sigGoodD0.thetaerr<90");
  pTNErrR->Draw();

  if (2 == nCharges) {
    can->cd(3);
    TProfile * pRPErrR = new TProfile("pRPErrR", "RP err vs. RPgen", 50, LOR, HIR, "s");
    tree.Project("pRPErrR", "dalitzHolderP.sigGoodD0.rhoerr:rhoPgen", "dalitzHolderP.sigGoodD0.rhoerr<0.4");
    pRPErrR->Draw();
    
    can->cd(4);
    TProfile * pTPErrR = new TProfile("pTPErrR", "TP err vs. RPgen", 50, LOR, HIR, "s");
    tree.Project("pTPErrR", "dalitzHolderP.sigGoodD0.thetaerr:rhoPgen", "dalitzHolderP.sigGoodD0.thetaerr<90");
    pTPErrR->Draw();
  }

  can->SaveAs("profilesVsRhoWithSpread.eps");

  // profile error dependences on rho showing error on mean:
  can = new TCanvas("errProfRhoE", "Errors vs. true rho, bars=error on mean", 800, 400 * nCharges); 
  can->Divide(2, nCharges);

  can->cd(1);
  TProfile * pRNErrRE = new TProfile("pRNErrRE", "RN err vs. RNgen", 50, LOR, HIR);
  tree.Project("pRNErrRE", "dalitzHolderN.sigGoodD0.rhoerr:rhoNgen", "dalitzHolderN.sigGoodD0.rhoerr<0.4");
  pRNErrRE->Fit("pol3");
  pRNErrRE->Draw();

  can->cd(2);
  TProfile * pTNErrRE = new TProfile("pTNErrRE", "TN err vs. RNgen", 50, LOR, HIR);
  tree.Project("pTNErrRE", "dalitzHolderN.sigGoodD0.thetaerr:rhoNgen", "dalitzHolderN.sigGoodD0.thetaerr<90");
  pTNErrRE->Fit("pol3");
  pTNErrRE->Draw();

  if (2 == nCharges) {
    can->cd(3);
    TProfile * pRPErrRE = new TProfile("pRPErrRE", "RP err vs. RPgen", 50, LOR, HIR);
    tree.Project("pRPErrRE", "dalitzHolderP.sigGoodD0.rhoerr:rhoPgen", "dalitzHolderP.sigGoodD0.rhoerr<0.4");
    pRPErrRE->Fit("pol3");
    pRPErrRE->Draw();
    
    can->cd(4);
    TProfile * pTPErrRE = new TProfile("pTPErrRE", "TP err vs. RPgen", 50, LOR, HIR);
    tree.Project("pTPErrRE", "dalitzHolderP.sigGoodD0.thetaerr:rhoPgen", "dalitzHolderP.sigGoodD0.thetaerr<90");
    pTPErrRE->Fit("pol3");
    pTPErrRE->Draw();
  }

  can->SaveAs("profilesVsRhoWithMeanError.eps");

  // profile error dependences on theta:
  can = new TCanvas("errDepTheta", "errDepTheta", 800, 400 * nCharges); 
  can->Divide(2, nCharges);

  can->cd(1);
  tree.Draw("dalitzHolderN.sigGoodD0.rhoerr:thetaNgen", "dalitzHolderN.sigGoodD0.rhoerr<0.4", "profs");

  can->cd(2);
  tree.Draw("dalitzHolderN.sigGoodD0.thetaerr:thetaNgen", "dalitzHolderN.sigGoodD0.thetaerr<90", "profs");

  if (2 == nCharges) {
    can->cd(3);
    tree.Draw("dalitzHolderP.sigGoodD0.rhoerr:thetaPgen", "dalitzHolderP.sigGoodD0.rhoerr<0.4", "profs");
    
    can->cd(4);
    tree.Draw("dalitzHolderP.sigGoodD0.thetaerr:thetaPgen", "dalitzHolderP.sigGoodD0.thetaerr<90", "profs");
  }

  can->SaveAs("profilesVsThetaWithSpread.eps");

  // Plot pulls:
  can = new TCanvas("pulls", "pulls", 800, 400 * nCharges); 
  can->Divide(2, nCharges);

  can->cd(1);
  TH1F * hpullRN = new TH1F("hpullRN", "hpullRN", 100, -7, 7);
  tree.Project("hpullRN", "(dalitzHolderN.sigGoodD0.rho-rhoNgen)/dalitzHolderN.sigGoodD0.rhoerr");
  hpullRN->Fit("gaus");
  hpullRN->Draw();

  can->cd(2);
  TH1F * hpullTN = new TH1F("hpullTN", "hpullTN", 100, -7, 7);
  tree.Project("hpullTN", "(dalitzHolderN.sigGoodD0.theta-thetaNgen)/dalitzHolderN.sigGoodD0.thetaerr");
  hpullTN->Fit("gaus");
  hpullTN->Draw();

  if (2 == nCharges) {
    can->cd(3);
    TH1F * hpullRP = new TH1F("hpullRP", "hpullRP", 100, -7, 7);
    tree.Project("hpullRP", "(dalitzHolderP.sigGoodD0.rho-rhoPgen)/dalitzHolderP.sigGoodD0.rhoerr");
    hpullRP->Fit("gaus");
    hpullRP->Draw();
    
    can->cd(4);
    TH1F * hpullTP = new TH1F("hpullTP", "hpullTP", 100, -7, 7);
    tree.Project("hpullTP", "(dalitzHolderP.sigGoodD0.theta-thetaPgen)/dalitzHolderP.sigGoodD0.thetaerr");
    hpullTP->Fit("gaus");
    hpullTP->Draw();
  }

  // plot the pulls in true rhoGen bins. For some reason, I'm not able to plot this directly 
  // (something happens to the canvas), so writing to a file and then reading back to plot:
  TFile * pullVsBinFile = new TFile(pullVsBinFileName, "RECREATE");

  TH1F * pullRNMean = new TH1F("pullRNMean", "RN pull mean vs RN", nBinsPulls, LOR, HIR);
  TH1F * pullRNSigma = new TH1F("pullRNSigma", "RN pull sigma vs RN", nBinsPulls, LOR, HIR);
  pullVsBin(pullRNMean, pullRNSigma, tree, "(dalitzHolderN.sigGoodD0.rho-rhoNgen)/dalitzHolderN.sigGoodD0.rhoerr", "rhoNgen");
  pullRNMean->Write();
  pullRNSigma->Write();

  TH1F * pullTNMean = new TH1F("pullTNMean", "TN pull mean vs RN", nBinsPulls, LOR, HIR);
  TH1F * pullTNSigma = new TH1F("pullTNSigma", "TN pull sigma vs RN", nBinsPulls, LOR, HIR);
  pullVsBin(pullTNMean, pullTNSigma, tree, "(dalitzHolderN.sigGoodD0.theta-thetaNgen)/dalitzHolderN.sigGoodD0.thetaerr", "rhoNgen");
  pullTNMean->Write();
  pullTNSigma->Write();

  if (2 == nCharges) {
    TH1F * pullRPMean = new TH1F("pullRPMean", "RP pull mean vs RP", nBinsPulls, LOR, HIR);
    TH1F * pullRPSigma = new TH1F("pullRPSigma", "RP pull sigma vs RP", nBinsPulls, LOR, HIR);
    pullVsBin(pullRPMean, pullRPSigma, tree, "(dalitzHolderP.sigGoodD0.rho-rhoPgen)/dalitzHolderP.sigGoodD0.rhoerr", "rhoPgen");
    pullRPMean->Write();
    pullRPSigma->Write();
    
    TH1F * pullTPMean = new TH1F("pullTPMean", "TP pull mean vs RP", nBinsPulls, LOR, HIR);
    TH1F * pullTPSigma = new TH1F("pullTPSigma", "TP pull sigma vs RP", nBinsPulls, LOR, HIR);
    pullVsBin(pullTPMean, pullTPSigma, tree, "(dalitzHolderP.sigGoodD0.theta-thetaPgen)/dalitzHolderP.sigGoodD0.thetaerr", "rhoPgen");
    pullTPMean->Write();
    pullTPSigma->Write();
  }

  pullVsBinFile->Close();
  delete pullVsBinFile;

  // read back and plot them:
  pullVsBinFile = new TFile(pullVsBinFileName);

  // first plot the pull means:
  can = new TCanvas("pullMeansVsRho", "pull means vs. rho", 
		    800, 400 * nCharges); 
  can->Divide(2 * nCharges, 1);

  can->cd(1);
  pullRNMean = (TH1F*)(pullVsBinFile->Get("pullRNMean"));
  pullRNMean->Fit("pol3");
  pullRNMean->Draw();

  can->cd(2);
  pullTNMean = (TH1F*)(pullVsBinFile->Get("pullTNMean"));
  pullTNMean->Fit("pol3");
  pullTNMean->Draw();

  if (2 == nCharges) {
    can->cd(3);
    pullRPMean = (TH1F*)(pullVsBinFile->Get("pullRPMean"));
    pullRPMean->Fit("pol3");
    pullRPMean->Draw();
    
    can->cd(4);
    pullTPMean = (TH1F*)(pullVsBinFile->Get("pullTPMean"));
    pullTPMean->Fit("pol3");
    pullTPMean->Draw();
  }

  can->SaveAs("pullMeansVsRho.eps");

  // now plot the pull sigmas:
  can = new TCanvas("pullSigmasVsRho", "pull sigmas vs. rho", 
		    800, 400 * nCharges); 
  can->Divide(2 * nCharges, 1);

  can->cd(1);
  pullRNSigma = (TH1F*)(pullVsBinFile->Get("pullRNSigma"));
  pullRNSigma->Fit("pol3");
  pullRNSigma->Draw();

  can->cd(2);
  pullTNSigma = (TH1F*)(pullVsBinFile->Get("pullTNSigma"));
  pullTNSigma->Fit("pol3");
  pullTNSigma->Draw();

  if (2 == nCharges) {
    can->cd(3);
    pullRPSigma = (TH1F*)(pullVsBinFile->Get("pullRPSigma"));
    pullRPSigma->Fit("pol3");
    pullRPSigma->Draw();
    
    can->cd(4);
    pullTPSigma = (TH1F*)(pullVsBinFile->Get("pullTPSigma"));
    pullTPSigma->Fit("pol3");
    pullTPSigma->Draw();
  }

  can->SaveAs("pullSigmasVsRho.eps");




  // Show the Gaussian fits from which these dependences were obtained:
  showPullVsBinPlots(nCharges);


  // Make plots of the dependences of the spreads of the errors on rho:
  can = new TCanvas("errs", "errs", 800, 400 * nCharges); 
  can->Divide(2, nCharges);

  can->cd(1);
  TH1F * hErrSpreadRN = plotSpread(pRNErrR, "hErrSpreadRN");
  hErrSpreadRN->Draw();
  
  can->cd(2);
  TH1F * hErrSpreadTN = plotSpread(pTNErrR, "hErrSpreadTN");
  hErrSpreadTN->Draw();
  
  if (2 == nCharges) {
    can->cd(3);
    TH1F * hErrSpreadRP = plotSpread(pRPErrR, "hErrSpreadRP");
    hErrSpreadRP->Draw();
    
    can->cd(4);
    TH1F * hErrSpreadTP = plotSpread(pTPErrR, "hErrSpreadTP");
    hErrSpreadTP->Draw();  
  }

  can->SaveAs("errorSpreadsVsRho.eps");

  // Make plots of the error distribution per ben:
  plotErrDistVsBin(tree, "dalitzHolderN.sigGoodD0.rhoerr", "rhoNgen", LOR, HIR);
  plotErrDistVsBin(tree, "dalitzHolderN.sigGoodD0.thetaerr", "rhoNgen", LOR, HIR);
  if (2 == nCharges) {
    plotErrDistVsBin(tree, "dalitzHolderP.sigGoodD0.rhoerr", "rhoPgen", LOR, HIR);
    plotErrDistVsBin(tree, "dalitzHolderP.sigGoodD0.thetaerr", "rhoPgen", LOR, HIR);
  }
}



//-----------------------------------------------------------
// Make and poly-fit the spread of a profile histogram of an error vs. rho.
// This will serve to tell us by how much to increase the error of a toy 
// experiment relative to the average for that given true-rho:
TH1F * plotSpread(const TProfile * prof, const char * name) {
  const int nBins = prof->GetNbinsX();

  TH1F * result = 
    new TH1F (name, 
	      TString("Spread of ") + TString(prof->GetTitle()) + TString(" vs. rho"), 
	      nBins, 
	      prof->GetBinLowEdge(1), 
	      prof->GetBinLowEdge(nBins+1)); // only way to get hist upper edge

  for (int b = 1; b <= nBins; ++b) {
    double sigma = prof->GetBinError(b);  // the spread was stored in the profile error
    double sigmaRMS = sigma / sqrt(2. * prof->GetBinEntries(b));  // formula for Gaussian sigma RMS

    result->SetBinContent(b, sigma);
    result->SetBinError(b, sigmaRMS);
  }
  
  result->Fit("pol3");
  return result;
}


//--------------------------------------------------------
// builds a binned name for a hist:
TString histName(const char * rawName, int bin) {
  TString name(rawName);
  name += "HistForBin";
  name += (float)bin;
  return name;
}


//--------------------------------------------------
pullVsBin(TH1F * hMean, TH1F * hSigma,
	  const TTree & tree, const char * var, const char * binVar) {
  // Make plots in of var in bins of binVar within LO and HI
  const int nBins = hMean->GetNbinsX();
  const double LO = hMean->GetBinLowEdge(1);
  const double HI = hMean->GetBinLowEdge(nBins+1);

  const double binWidth = (HI-LO)/nBins;

  double mean = 0;
  double meanErr = 0;
  double sigma = 0;
  double sigmaErr = 0;

  for (int b = 0; b < nBins; ++b) {
    double min = LO + b * binWidth;
    double max = HI + b * binWidth;

    TH1F hist(histName(hMean->GetName(), b), histName(hMean->GetName(), b), 70, -7, 7);

    TString cut = binCut(binVar, min, max);

    tree.Project(hist.GetName(), var, cut);
    hist.Fit("gaus");
    hist.Write();

    mean = hist.GetFunction("gaus")->GetParameter(1);
    meanErr = hist.GetFunction("gaus")->GetParError(1);

    sigma = hist.GetFunction("gaus")->GetParameter(2);
    sigmaErr = hist.GetFunction("gaus")->GetParError(2);

    hMean->SetBinContent(b+1, mean);

    hMean->SetBinError(b+1, meanErr);

    hSigma->SetBinContent(b+1, sigma);
    hSigma->SetBinError(b+1, sigmaErr);
  }
}


TString binCut(const char * binVar, double min, double max) {
  // Make a bin cut:
  TString cut(binVar);
  cut += "<";
  cut += max;
  cut += "&&";
  cut += binVar;
  cut += ">";
  cut += min;
  return cut;
}


void plotErrDistVsBin(const TTree & tree, const char * var, const char * binVar, double LO, double HI) {
  // Plot the fit error distribution in each rho bin:
  int nX = 5;
  int nY = nBinsPulls / nX;
  if (nY * nX < nBinsPulls) {
    ++nY;
  }

  const double binWidth = (HI-LO)/nBinsPulls;
  
  TString canName = TString(var) + TString("-errDistPerBin");

  TCanvas * can = new TCanvas(canName, canName, 200 * nX, 200 * nY);
  can->Divide(nX, nY);
  for (int b = 0; b < nBinsPulls; ++b) {
    double min = LO + b * binWidth;
    double max = HI + b * binWidth;

    TString cut = binCut(binVar, min, max) + +TString("&&") + TString(var) + TString("<90");
    can->cd(b+1);
    tree.Draw(var, cut);
  }
  can->SaveAs(TString(var) + "-error-dist-in-each-bin.eps");
}


// Show the pull distributions on each bin for the pull plots created by pullVsBin(...):
void showPullVsBinPlots(int nCharges) {
  TFile * pullVsBinFile = new TFile(pullVsBinFileName);
  const int nX = 5;
  const int nY = 4;
  if (nBinsPulls != nX * nY) {
    cerr << "********************************************" << endl;
    cerr << "**** Note: Not all bins are being shown ****" << endl;
    cerr << "********************************************" << endl;
  }

  getPlotPulls(pullVsBinFile, "pullRNMean", nX, nY);
  getPlotPulls(pullVsBinFile, "pullTNMean", nX, nY);
  if (2 == nCharges) {
    getPlotPulls(pullVsBinFile, "pullRPMean", nX, nY);
    getPlotPulls(pullVsBinFile, "pullTPMean", nX, nY);
  }
}

void getPlotPulls(const TFile * file, const char * rawName, int nX, int nY) {
  TCanvas * can = new TCanvas(rawName, rawName, 200 * nX, 200 * nY);
  can->Divide(nX, nY);
  for (int b = 0; b < nBinsPulls; ++b) {
    TH1F * hist = (TH1F*)(file->Get(histName(rawName, b)));
    can->cd(b+1);
    hist->Draw();
  }

  can->SaveAs(TString(rawName) + TString("-pullDistPerBin.eps"));
}

  
