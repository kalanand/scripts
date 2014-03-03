// Functions in this file:
//
// combineFiles(h3writeDir + h3File) is used to combine a TH3F from
// all the files produced by CLTreeMaker::makeCLHist3() into a single
// file.
//
// makeProjs(h3readDir + h3File) reads that TH3F from the single file
// and makes 1d and 2d projections, then writes them to a local root
// file.
//
// readProjHists(h3readDir + projFileName) reads the histograms back
// from the file.
//
// plotProjs() plots them
//
//
// Functions that call the above functions for all configurations
// ("All", "Stat", "StatSyst", "Bias") and take 4 corresponding Bool_t's:
// combineFilesAll(Bool_t all, Bool_t stat, Bool_t statSyst, Bool_t bias)
// makeProjsAll(Bool_t all, Bool_t stat, Bool_t statSyst, Bool_t bias)
// plotProjsAll(Bool_t all, Bool_t stat, Bool_t statSyst, Bool_t bias)
//

const TString h3readDir = "/nfs/farm/babar/AWG16/Breco/abi/";
const TString h3writeDir = "/afs/slac.stanford.edu/g/babar/work/a/abi/";

const TString h3File = "makeCLHist3-all.root";
const TString h3FileBias = "makeCLHist3-bias-all.root";
const TString h3FileStat = "makeCLHist3-stat-all.root";
const TString h3FileStatSyst = "makeCLHist3-statsyst-all.root";

const TString projFileName = "makeCLHists-projections.root";
const TString projFileNameBias = "makeCLHists-projections-bias.root";
const TString projFileNameStat = "makeCLHists-projections-stat.root";
const TString projFileNameStatSyst = "makeCLHists-projections-statsyst.root";

const TString h3HistName = "hCLAll-sum";
const TString h3HistNameBias = "hCLAllBias-sum";
const TString h3HistNameStat = "hCLStat-sum";
const TString h3HistNameStatSyst = "hCLStatSyst-sum";

TFile * projFile = 0;

CLHistProjector * proj = 0;

TH1F * hR = 0;
TH1F * hG = 0;
TH1F * hD = 0;

TH2F * hRG = 0;
TH2F * hRD = 0;
TH2F * hGD = 0;


const char * rawFiles[20] = 
{"/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-1.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-2.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-3.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-4.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-5.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-6.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-7.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-8.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-9.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-10.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-11.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-12.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-13.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-14.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-15.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-16.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-17.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-18.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-19.root",
 "/nfs/farm/babar/AWG16/Breco/abi/makeCLHist3-20.root"
};



// Need to give first argument by hand
void makeProjs(const char * fileName, // set to   h3readDir + h3File
	       const TString fileToSaveHists = projFileName,
	       const TString histName = h3HistName) {
  
  cout << "reading hist3" << endl;
  TFile * file = new TFile(fileName);
  TH3F * clAll = (TH3F*)file->Get(histName);

  if (0 == clAll) {
    cout << "null histogram read. Quitting." << endl;
    return;
  }

  cout << "making projections" << endl;
  proj = new CLHistProjector(clAll);
  
  proj->makeAllProj();

  cout << "Saving projections" << endl;
  proj->saveHists(fileToSaveHists);
}


void makeProjsAll(Bool_t doAll, 
		  Bool_t doStat, 
		  Bool_t doStatSyst, 
		  Bool_t doBias) {
  if (doAll){
    makeProjs(h3readDir + h3File,
	      h3writeDir + projFileName,
	      h3HistName);
  }
  
  if (doStat){
    makeProjs(h3readDir + h3FileStat,
	      h3writeDir + projFileNameStat,
	      h3HistNameStat);
  }
  
  if (doStatSyst){
    makeProjs(h3readDir + h3FileStatSyst,
	      h3writeDir + projFileNameStatSyst,
	      h3HistNameStatSyst);
  }
  
  if (doBias){
    makeProjs(h3readDir + h3FileBias,
	      h3writeDir + projFileNameBias,
	      h3HistNameBias);
  }
}




void format(TH1 * hist, const char * x, const char * y) {
  hist->SetXTitle(x);
  hist->SetYTitle(y);
  hist->UseCurrentStyle();

  hist->GetXaxis()->SetTitleOffset(0.01);
  hist->SetLabelOffset(0.02.,"X");
  hist->SetLabelOffset(0.02.,"Y");

  hist->SetNdivisions(-504,"X");
  hist->SetNdivisions(-506,"Y");
}


void plotProjs(const char * fileName = 0,
	       const char * name = 0) { // if !0, first reads the file
  if (0 != fileName && 0 != name) {
    readProjHists(fileName, name);
  }

  cout << "making plots" << endl;
  
  gROOT->SetStyle("BABAR");
  int colors[3] = {kBlue, kGreen, kYellow};
  gStyle->SetPalette(3, colors);

  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.20);

  gStyle->SetHistMinimumZero();
  
  gStyle->SetTitleXSize(0.1);
  gStyle->SetTitleYSize(0.1);

  gStyle->SetLabelSize(0.08, "X");
  gStyle->SetLabelSize(0.08, "Y");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(10);


  format(hR, "r_{B}",    "1-#alpha");
  format(hG, "#gamma (deg.)", "1-#alpha");
  format(hD, "#delta (deg.)", "1-#alpha");

  format(hRG, "r_{B}",    "#gamma (deg.)");
  format(hRD, "r_{B}",    "#delta (deg.)");
  format(hGD, "#gamma (deg.)", "#delta (deg.)");

  // Some changes for the 1D plot formats:  
  hR->SetNdivisions(506,"Y");
  hG->SetNdivisions(506,"Y");
  hD->SetNdivisions(506,"Y");


  const int NLEVELS = 3;
  int colors[NLEVELS] = {kBlue, kGreen, kYellow};
  double levels[NLEVELS] = {0.0, 0.199, 0.739};

  // So as not to get an ugly legend, need to have the third level
  // of hGD be the highest bin. So find the highest bin:
  double maxCLGD = 0;
  for (int iG = 1; iG < hGD->GetNbinsX(); ++iG) {
    for (int iD = 1; iD < hGD->GetNbinsY(); ++iD) {
      if (maxCLGD < hGD->GetBinContent(iG, iD)) {
	maxCLGD = hGD->GetBinContent(iG, iD);
      }
    }
  }
   
  double levelsGD[NLEVELS] = {0.0, 0.199, maxCLGD+1.e-08};

  gStyle->SetPalette(NLEVELS, colors);
  
  hRG->SetContour(NLEVELS, levels);
  hRD->SetContour(NLEVELS, levels);
  hGD->SetContour(NLEVELS, levelsGD);
  
  TCanvas * can = new TCanvas("c", "c", 1200, 800);
  can->Divide(3,2);
  
  // The 2D plots need more space on the right for the legend:
  gStyle->SetPadRightMargin(0.1);
  can->cd(1);
  hRG->Draw("colz");
  can->cd(2);
  hRD->Draw("colz");
  can->cd(3);
  hGD->Draw("colz");
  gStyle->SetPadRightMargin(0.05);
  can->cd(4);
  hR->Draw();
  can->cd(5);
  hG->Draw();
  can->cd(6);
  hD->Draw();

  cout << "Saving eps file" << endl;
  TString saveName = name;
  can->SaveAs(saveName + "-CL-proj.eps");

  // Make the 2D plot for the paper:
  TCanvas * can2D = new TCanvas("can2D", "can2D", 400, 400);
  can2D->GetPad(0)->SetRightMargin(0.15);
  hRG->Draw("colz");
  TPaletteAxis * pal = (TPaletteAxis*)hRG->FindObject("palette");
  // Adjust the position and label size of the palette axis ("legend"):
  pal->SetLabelSize(0.06);
  pal->SetX1NDC(0.9/pal->GetX1() * pal->GetX1NDC());
  pal->SetX2NDC(0.9/pal->GetX1() * pal->GetX2NDC());
  can2D->SaveAs(saveName + "-CL-proj-RG.eps");


  // Make the 1D plot for the paper:
  TCanvas * can1D = new TCanvas("can1D", "can1D", 400, 400);
  hG->Draw();
  can1D->SaveAs(saveName + "-CL-proj-G.eps");
}


// Run plotProjs on all:
plotProjsAll(Bool_t doAll, Bool_t doStat, Bool_t doStatSyst, Bool_t doBias) {
  if (kTRUE == doAll){
    plotProjs(h3readDir + projFileName, h3HistName);
  }

  if (kTRUE == doStat){
    plotProjs(h3readDir + projFileNameStat, h3HistNameStat);
  }

  if (kTRUE == doStatSyst){
    plotProjs(h3readDir + projFileNameStatSyst, h3HistNameStatSyst);
  }

  if (kTRUE == doBias){
    plotProjs(h3readDir + projFileNameBias, h3HistNameBias);
  }
}



// Sum all the histograms in a directory. Here, need to give the 
// first argument explicitly, otherwise there is a crash due to 
// TString and how the interpreter handles temporary variables,
// or something like that:
void combineFiles(const char * writeFile,  // set to   h3writeDir + h3File,
		  const char * histName = "hCLAll",
		  const char * files[] = rawFiles, 
		  int nFiles = 20) {   

  proj = new CLHistProjector();
  proj->combine(files, nFiles, histName);

  TFile * file = new TFile(writeFile, "RECREATE");
  if (0 == file || (0 != file && kFALSE == file->IsOpen())) {
    cout << "Can't open file " << writeFile << endl;
  }

  cout << "Saving the new h3 to " << writeFile << endl;
  proj->h3()->Write();
  file->Close();
  delete file;
}



void combineFilesAll(Bool_t doAll, 
		     Bool_t doStat, 
		     Bool_t doStatSyst, 
		     Bool_t doBias) {
  if (doAll){
    combineFiles(h3readDir + h3File, "hCLAll");
  }
  
  if (doStat){
    combineFiles(h3readDir + h3FileStat, "hCLStat");
  }
  
  if (doStatSyst){
    combineFiles(h3readDir + h3FileStatSyst, "hCLStatSyst");
  }
  
  if (doBias){
    combineFiles(h3readDir + h3FileBias, "hCLAllBias");
  }
}





void readProjHists(const char * fileName, // set to h3readDir + projFileName,
		   const char * name = h3HistName) {
  cout << "readProjHists called to read " << fileName 
       << endl;
  
  projFile = new TFile(fileName);

  hR = (TH1F*)projFile->Get(TString(name) + TString("-rB"));
  hG = (TH1F*)projFile->Get(TString(name) + TString("-gamma"));
  hD = (TH1F*)projFile->Get(TString(name) + TString("-delta"));
  
  hRG = (TH2F*)projFile->Get(TString(name) + TString("-rB-gamma"));
  hRD = (TH2F*)projFile->Get(TString(name) + TString("-rB-delta"));
  hGD = (TH2F*)projFile->Get(TString(name) + TString("-gamma-delta"));
}
  


