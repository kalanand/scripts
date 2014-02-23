
#include <iostream>
using namespace std;

#include <assert.h>

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

#include "BToDKTo3piK/CLHistProjector.hh"

ClassImp(CLHistProjector)

CLHistProjector::CLHistProjector(const TH3F * h3) :
  _h3(h3),
  _ownH3(false)
{
}
 

CLHistProjector::~CLHistProjector() {
  delete _hRG;
  delete _hRD;
  delete _hGD;
  delete _hR;
  delete _hG;
  delete _hD;
  deleteH3();
}


void CLHistProjector::deleteH3() {
  if (_ownH3) {
    TH3F * h3Own = (TH3F*)_h3;
    delete h3Own;
  }
}



const char * CLHistProjector::axisName(const TAxis * axis) const {
  if (axisR() == axis) {
    return nameR();
  }
  if (axisG() == axis) {
    return nameG();
  }
  if (axisD() == axis) {
    return nameD();
  }
  return 0;
}


// Translate the i# bin numbers into iR, iG, iD bin numbers, given the axes:
void CLHistProjector::translate(const TAxis * axes[3], 
				int indices[3],
				int & iR, int & iG, int & iD) const {
  for (int j = 0; j < 3; ++j) {
    if (axisR() == axes[j]) {iR = indices[j];}
    if (axisG() == axes[j]) {iG = indices[j];}
    if (axisD() == axes[j]) {iD = indices[j];}
  }
}



// "Project" the 3D hist onto 1D: 
TH1F * CLHistProjector::projOnto(const TAxis * axis) const { 
  // make the resulting histogram:
  TString name = _h3->GetName() + TString("-") + axisName(axis);

  int nBins = axis->GetNbins();
  double lo = axis->GetXmin();
  double hi = axis->GetXmax();

  TH1F * result = new TH1F(name, name, nBins, lo, hi);

  // find the other axis:
  const TAxis * axis2 = 0;
  const TAxis * axis3 = 0;

  if (axisR() == axis) {
    axis2 = axisG();
    axis3 = axisD();
  }
  else if (axisG() == axis) {
    axis2 = axisD();
    axis3 = axisR();
  }
  else if (axisD() == axis) {
    axis2 = axisR();
    axis3 = axisG();
  }

  // Prepare for use later:
  const TAxis * axes[3] = {axis, axis2, axis3};

  // For each bin in the result histogram:
  for (int i1 = 1; i1 <= axis->GetNbins(); ++i1) {
    result->SetBinContent(i1, 2.0);  // start from ~2 * max value
    // loop over bins of the other axes:
    for (int i2 = 1; i2 <= axis2->GetNbins(); ++i2) {
      for (int i3 = 1; i3 <= axis3->GetNbins(); ++i3) {
	// Translate the i# bin numbers into iR, iG, iD bin numbers:
	int indices[3] = {i1, i2, i3};
	int iR, iG, iD;
	translate(axes, indices, iR, iG, iD);

	// and keep the smallest _h3-valued bin:
	double h3Value = _h3->GetBinContent(iR, iG, iD);
	if (result->GetBinContent(i1) > h3Value) {
	  result->SetBinContent(i1, h3Value);
	}
      }
    }
  }

  return result;
}


TH2F * CLHistProjector::projOnto(const TAxis * axis1, 
				 const TAxis * axis2) const {
  // make the resulting histogram:
  TString name = _h3->GetName() + TString("-") + axisName(axis1)
    + TString("-") + axisName(axis2);

  int nBins[2] = {axis1->GetNbins(), axis2->GetNbins()};
  double lo[2] = {axis1->GetXmin(), axis2->GetXmin()};
  double hi[2] = {axis1->GetXmax(), axis2->GetXmax()};

  TH2F * result = new TH2F(name, name, 
			   nBins[0], lo[0], hi[0],
			   nBins[1], lo[1], hi[1]);

  // find the other axis:
  const TAxis * axis3 = axisR();

  if (axisR() == axis1 || axisR() == axis2) {
    axis3 = axisG();
    if (axisG() == axis1 || axisG() == axis2) {
      axis3 = axisD();
    }
  }

  // Prepare for use later:
  const TAxis * axes[3] = {axis1, axis2, axis3};

  // For each bin in the result histogram:
  for (int i1 = 1; i1 <= axis1->GetNbins(); ++i1) {
    for (int i2 = 1; i2 <= axis2->GetNbins(); ++i2) {
      result->SetBinContent(i1, i2, 2.0);  // start from ~2 * max value
      // loop over bins of the other axis:
      for (int i3 = 1; i3 <= axis3->GetNbins(); ++i3) {
	// Translate the i# bin numbers into iR, iG, iD bin numbers:
	int indices[3] = {i1, i2, i3};
	int iR, iG, iD;
	translate(axes, indices, iR, iG, iD);

	// and keep the smallest _h3-valued bin:
	double h3Value = _h3->GetBinContent(iR, iG, iD);
	if (result->GetBinContent(i1, i2) > h3Value) {
	  result->SetBinContent(i1, i2, h3Value);
	}
      }
    }
  }

  return result;
}


void CLHistProjector::makeAllProj() {
  _hR = projOnto(axisR());
  _hG = projOnto(axisG());
  _hD = projOnto(axisD());

  _hRG = projOnto(axisR(), axisG());
  _hRD = projOnto(axisR(), axisD());
  _hGD = projOnto(axisG(), axisD());
}



void CLHistProjector::combine(const char * files[], 
			      int nFiles,
			      const char * histName) {
  // Combines all the histograms named histName in the files into a single
  // histogram and puts it in _h3.
  // Note that the histogram for each file has a different range of 
  // rB bins, which is taken into account:

  int nBinsRPerHist = 0;
  int nBinsR = 0;
  int nBinsG = 0;
  int nBinsD = 0;

  double loR = 10;  // the R limits have "inverse" initial values 
  double hiR = 0;
  double loG = -1;  // the other limits have something arbitrary
  double hiG = -1;
  double loD = -1;
  double hiD = -1;

  // Read the files once to determine their limits and nbins.
  // The algorithm closes unused files, so as not to hold too much
  // stuff in memory when we create the big _h3.
  cout << "Reading files to find limits:" << endl;
  int f;
  for (f = 0; f < nFiles; ++f) {
    TFile file(files[f]);   // not a pointer, so lives only in this loop
    const TH3F * h3 = (const TH3F*)file.Get(histName);

    if (0 == h3) {
      cerr << "File " << files[f] << " contains no object named "
	   << histName  << ". It does include: " << endl;
      file.ls();
      cerr << "Skipping this file." << endl;
      continue;
    }

    const TAxis * rAxis = h3->GetXaxis();
    const TAxis * gAxis = h3->GetYaxis();
    const TAxis * dAxis = h3->GetZaxis();

    cout << " " << files[f] << ", rB limits= " 
	 << rAxis->GetXmin() << " -- " << rAxis->GetXmax() << endl;

    // Find the rB limits from the extremes of all histograms:
    if (loR > rAxis->GetXmin()) {loR = rAxis->GetXmin();}
    if (hiR < rAxis->GetXmax()) {hiR = rAxis->GetXmax();}

    if (0 == f) {
      // on the first file, calculate the number of rB bins:
      nBinsRPerHist = rAxis->GetNbins();
      nBinsR = nBinsRPerHist * nFiles;
      nBinsG = gAxis->GetNbins();
      nBinsD = dAxis->GetNbins();

      // Also get the delta and gamma histogram limits: 
      loG = gAxis->GetXmin();
      hiG = gAxis->GetXmax();

      loD = dAxis->GetXmin();
      hiD = dAxis->GetXmax();
    }
    else {
      // on subsequent files, just check that the # of bins hasn't changed:
      assert(rAxis->GetNbins() == nBinsRPerHist);
      assert(gAxis->GetNbins() == nBinsG);
      assert(dAxis->GetNbins() == nBinsD);

      // For gamma and delta, also the limits shouldn't change. Not the same
      // for rB, whose range is determined from the extremes of all 
      // histograms:
      assert(gAxis->GetXmin() == loG);
      assert(dAxis->GetXmin() == loD);

      assert(gAxis->GetXmax() == hiG);
      assert(dAxis->GetXmax() == hiD);
    }
  }

  cout << "Found nBins: R=" << nBinsR << ", G=" << nBinsG << ", D=" << nBinsD
       << endl
       << "  low Range:  R=" << loR << ", G=" << loG << ", D=" << loD
       << endl
       << " high Range:  R=" << hiR << ", G=" << hiG << ", D=" << hiD
       << endl;
  
  // Done getting limits. Create the resulting histogram:
  deleteH3();

  TString h3Name = histName;
  h3Name += "-sum";

  TH3F * h3new = new TH3F(h3Name, h3Name, 
			  nBinsR, loR, hiR,
			  nBinsG, loG, hiG,
			  nBinsD, loD, hiD);
    
  // Extract the histograms from the files and fill into h3new:
  for (f = 0; f < nFiles; ++f) {
    cout << "Adding histogram from " << files[f] << endl;

    TFile file(files[f]);   // not a pointer, so lives only in this loop
    const TH3F * h3toAdd = (const TH3F*)file.Get(histName);

    const TAxis * axisR = h3toAdd->GetXaxis();
    int nBinsR = axisR->GetNbins();

    const TAxis * axisG = h3toAdd->GetYaxis();
    int nBinsG = axisG->GetNbins();

    const TAxis * axisD = h3toAdd->GetZaxis();
    int nBinsD = axisD->GetNbins();

    // Loop over all bins in h3toAdd and fill their contents into h3new:
    for (int iR = 1; iR <= nBinsR; ++iR) {
      double rB = axisR->GetBinCenter(iR);
      for (int iG = 1; iG <= nBinsG; ++iG) {
	double gamma = axisG->GetBinCenter(iG);
	for (int iD = 1; iD <= nBinsD; ++iD) {
	  double delta = axisD->GetBinCenter(iD);

	  int filledBin = h3new->Fill(rB, gamma, delta, 
				      h3toAdd->GetBinContent(iR, iG, iD));

	  // Error checking that bin hasn't been filled twice:
	  if(h3new->GetBinContent(filledBin) 
	     != h3toAdd->GetBinContent(iR, iG, iD)) {
	    cout << "h3new's bin " << filledBin << " was already filled. iR="
		 << iR << ", iG=" << iG << "iD=" << iD << endl;
	    assert(0);
	  }
	}
      } 
    }  // end loop over 3D bins
  }  // end loop over files
  // Copy the mutable h3New into the const _h3:
  _h3 = h3new;
  _ownH3 = true;
}



void CLHistProjector::saveHists(const char * fileName, bool saveH3) const {
  TFile file(fileName, "RECREATE");
  _hRG->Write();
  _hRD->Write();
  _hGD->Write();
  _hR->Write();
  _hG->Write();
  _hD->Write();

  if (saveH3) {
    _h3->Write();
  }

  cout << "Saved histograms to file " << fileName << endl;
}
