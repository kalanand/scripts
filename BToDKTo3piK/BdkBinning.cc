#include <assert.h>
#include <ostream>
using namespace std;

#include "BToDKTo3piK/BdkDP.hh"
#include "BToDKTo3piK/BdkIsoOp.hh"
#include "BToDKTo3piK/BdkSymOp.hh"
#include "BToDKTo3piK/BdkBinWgtList.hh"
#include "BToDKTo3piK/BdkBinning.hh"

#include "TString.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"

ClassImp(BdkBinning)

BdkBinning::BdkBinning(int nDivMid,
		       TTree * dataSR,
		       TTree * dataSB) :
  _wgtListsMade(kFALSE),
  _wfCalculated(kFALSE),
  _evErr(TMatrixD(BdkIsoOp::NOPS, BdkIsoOp::NOPS))
{
  setNDivMid(nDivMid);
  setData(dataSR, dataSB);
}


BdkBinning::~BdkBinning() {}


void BdkBinning::setNDivMid(int nDivMid) {
  _nDivMid = nDivMid;

  // Remove any old bins:
  _bins.clear();

  // Make the bins. We need at least 2*_nDivMid bins, so make 3*:
  const int nDiv = 3 * _nDivMid;
  const double space = BdkDP::symValue() / _nDivMid;
  int lastNBins = 0;

  cout << "BdkBinning::setNDivMid(" << nDivMid << "): Creating bins" << endl;
  for (int ix = 0; ix <= nDiv; ++ix) {
    for (int iy = 0; iy <= nDiv; ++iy) {
      // At each grid point put the corner of a straight bin and a mirror bin:
      BdkP2 point(ix * space, iy * space);
      BdkTriBin bin    (point, space, kFALSE);
      BdkTriBin binMirr(point, space, kTRUE );
      addBin(bin);
      addBin(binMirr);
    } // end loop on column
    cout << "Added " << _bins.size() - lastNBins 
	 << " bins to vertical column " << ix << endl;
    
    lastNBins = _bins.size();
  } // end loop on row
 
  _wgtListsMade = kFALSE;
  _wfCalculated = kFALSE;
}


void BdkBinning::setData(TTree * dataSR, TTree * dataSB){
  _dataSR = dataSR;
  _dataSB = dataSB;
  _wfCalculated = kFALSE;
}


// Divide a bin into three. Returns true if found the bin, false otherwise:
Bool_t BdkBinning::divideBin(int index) {
  // First, make all the symmetry-related bins:
  BdkTriBin & bin = _bins.at(index);
  vector<BdkTriBin> symBins;
  symBins.push_back(bin);
  symBins.push_back(bin.R());
  symBins.push_back(bin.RR());
  symBins.push_back(bin.E());
  symBins.push_back(bin.ER());
  symBins.push_back(bin.ERR());

  for (int b = 0; b < symBins.size(); ++b) {
    // Find the corresponding bin in _bins:
    int oldBinIndex = binIndex(symBins[b]);
    BdkTriBin & oldBin = _bins.at(oldBinIndex);

    // Divide it:
    vector<BdkTriBin> newBins = oldBin.divide();

    // Remove oldBin and add its new sub-bins:
    vector<BdkTriBin>::iterator iter = _bins.begin();
    for (int dummy = 0; dummy < oldBinIndex; ++dummy) {
      ++iter;
    }
    _bins.erase(iter);

    for (int nb = 0; nb < newBins.size(); ++nb){
      addBin(newBins[nb]);
    }
  } // end loop on the symmetry-related bins of _bins[index]
    
  _wgtListsMade = kFALSE;
  _wfCalculated = kFALSE;
  
  return kTRUE;  // success
}


Bool_t BdkBinning::divideBin(const BdkP2 & point) {
  int index = binIndexContaining(point);
  if (0 > index) {
    return kTRUE;
  }
  return divideBin(index);
}


// Builds the symmetry relations between the bins:
void BdkBinning::buildWgtLists() {
  if (kTRUE == _wgtListsMade) {
    return;
  }

  cout << "BdkBinning::buildWgtLists(): making weight lists for " 
       <<  _bins.size() << " bins:" << endl;

  // Get each bin's symmetry-related bins:
  for(int b = 0; b < _bins.size(); ++b) {
    if (b%10 == 0) { // monitor progress
      cout << "Bin " << b << endl;
    }

    // The bin whose symmetry related bins we're looking for
    BdkTriBin & theBin = _bins[b];

    // Its symmetry-related bins (temp objects):
      const BdkTriBin & binR = _bins[binIndex(theBin.R())];
      const BdkTriBin & binRR = _bins[binIndex(theBin.RR())];
      const BdkTriBin & binE = _bins[binIndex(theBin.E())];
      const BdkTriBin & binER = _bins[binIndex(theBin.ER())];
      const BdkTriBin & binERR = _bins[binIndex(theBin.ERR())];

    // For each of the operator-related binWgtList objects of this bin,
    // set all the symmetry-related bins:
    for (int op = 0; op < BdkIsoOp::NOPS; ++op) {
      BdkBinWgtList & wgtList = theBin._wgtList[op];
      
      // Find the bins in _bins that are equivalent to theBin's
      // symmetry-related bins:
      wgtList.setI(theBin);      
      wgtList.setR(binR);      
      wgtList.setRR(binRR);      
      wgtList.setE(binE);      
      wgtList.setER(binER);      
      wgtList.setERR(binERR);
    }
  }
  cout << endl;
  _wgtListsMade = kTRUE;
}


void BdkBinning::setDataToBins(TTree * tree) {
  cout << "BdkBinning::setDataToBins(...): Filling bins from tree \"" 
       << tree->GetName() << "\"." << endl;

  if (tree != _dataSB && tree != _dataSR) {
    cout << "*** BdkBinning::setDataToBins(): given unfamiliar tree."
	 << endl;
    assert(0);
  }
  
  double s12, s13;
  TBranch * bS12 = tree->GetBranch("s12");
  TBranch * bS13 = tree->GetBranch("s13");
  bS12->SetAddress(&s12); 
  bS13->SetAddress(&s13); 
  
  int nRejectedPoints = 0;
  // loop over events:
  for (int e = 0; e < tree->GetEntries(); ++e) {
    bS12->GetEntry(e);
    bS13->GetEntry(e);
    BdkP2 point(s12, s13);
    
    if (kFALSE == point.inDalitz()) {
      ++nRejectedPoints;
      continue;  // don't take points outside the DP
    }
    
    // Loop over bins and find the one that contains the event:
    for (int b = 0; b < _bins.size(); ++b) {
      if (_bins[b].contains(point)) {
	if (tree == _dataSR) {
	  _bins[b].accumNEventsSR();
	}
	else {
	  _bins[b].accumNEventsSB();
	}
	break; // out of the loop on bins and on to the next event
      } // end if this is the bin containing the point
    } // end loop on bins
  } // end loop on events

  cout << tree->GetEntries() - nRejectedPoints
       << " points added to bins. " << nRejectedPoints
       << " points outside the DP were rejected." << endl;
}



void BdkBinning::calcWFs(double mDratioSRtoSB) {
  if (kTRUE == _wfCalculated) {
    return;
  }

  // Read the data sets:
  setDataToBins(_dataSR);
  setDataToBins(_dataSB);

  // then calculate the wf's:
  for (int b = 0; b < _bins.size(); ++b) {
    _bins[b].calcWaveFunc(mDratioSRtoSB);
  }

  _wfCalculated = kTRUE;
}


// Calculate operator expectation values and errors:
void BdkBinning::calcEVs() {
  // set up the necessary information:
  calcWFs(); 
  buildWgtLists();

  // calculate expectation values:
  int b;  // bin index

  // Calculate the sum over the squares of the wave functions, which
  // is the denominator in the EV expression:
  double sumWFSQ = 0;
  for (b = 0; b < _bins.size(); ++b) {
    double wf = _bins[b].waveFunc();
    sumWFSQ += wf * wf;
  }

  // Loop over operators and calcuate the EV's:
  for (int isoOp = 0; isoOp < BdkIsoOp::NOPS; ++isoOp) {
    // Zero-out the EV:
    _ev[isoOp] = 0;
    
    // Calculate the numerator:
    for (b = 0; b < _bins.size(); ++b) {
      double wf = _bins[b].waveFunc();
  
      // This wave function has to be multiplied by the sum of the
      // wave functions and weights of all the bins with a symmetry
      // relation to this bin:
      double sumWgtsWFs = 0;
      const BdkBinWgtList & wgtList = _bins[b].wgtList(isoOp);
      for (int symOp = 0; symOp < BdkSymOp::NOPS; ++symOp) {
	const BdkBinWgt & binWgt = wgtList.binWgt(symOp);
	double binWeight = binWgt.weight();
	double binWf = binWgt.bin()->waveFunc();
	sumWgtsWFs += binWeight * binWf;
      }
     
      // Do the said multiplication:
      _ev[isoOp] += wf * sumWgtsWFs;
    }

    // Divide by the denominator:
    _ev[isoOp] /= sumWFSQ;
  } // End loop on isoOp, calculation of expectation values

  // Now calculate the error matrix:
  int i;
  for (i = 0; i < BdkIsoOp::NOPS; ++i) { // loop on matrix rows
    for (int j = i; j < BdkIsoOp::NOPS; ++j) {  // on matrix columns
      _evErr[j][j] = 0.0; // initialize
      for (int a = 0; a < _bins.size(); ++a) {  // on bins
	const BdkTriBin & bin = _bins[a];
	const double wf = bin.waveFunc();
	const double wfErr = bin.waveFuncErr();

	// Get the two wgtLists and calculate the first terms (sums)
	// inside the square brackets of Eq (45) in isospin-detail-13:
	const BdkBinWgtList * wgtList[2] = {&(bin.wgtList(i)), 
					    &(bin.wgtList(j))};
					    
	double term[2] = {0,0};
	// Loop over the two lists to calculate the two terms:
	for (int f = 0; f < 2; ++f) {
	  // For each list, loop over the symmetry-related bins and
	  // get the weights:
	  for (int symOp = 0; symOp < BdkSymOp::NOPS; ++symOp) {
	    // Get the binWgt of bin with its symmetry-related bin:
	    const BdkBinWgt & binWgt = wgtList[f]->binWgt(symOp);

	    // From this, get the weight:
	    const double binWeight = binWgt.weight();

	    // and also get the symmetry related bin:
	    const BdkTriBin * symBin = binWgt.bin();

	    // Get the symmetry-related bin's wave function:
	    const double binWf = symBin->waveFunc();

	    // Get the reverse weight, i.e., the symmetry-related
	    // bin's weight for to transform into bin:
	    const BdkBinWgtList & reverseWgtList = 
	      symBin->wgtList(wgtList[f]->op());

	    const BdkBinWgt * reverseBinWgt = reverseWgtList.binWgt(bin);
	    assert(reverseBinWgt != 0);
	    const double reverseWgt = reverseBinWgt->weight();

	    // Now make the necessary calculation:
	    term[f] += binWf * (binWeight + reverseWgt);
	  } // end loop on symOp to calculate terms
	 
	} // end loop on f (the two square brackets in Eq (45)) 
	
	// Now subtract from term[f] the second term in Eq (45):
	term[0] -= 2 * wf * _ev[i];
	term[1] -= 2 * wf * _ev[j];

	// Complete the product of Eq (45) and accumulate the term
	// into element of matrix:
	_evErr[i][j] += term[0] * wfErr * wfErr * term[1] / sumWFSQ;
      } // end loop on bins (index a)
    } // end loop on matrix columns
  } // end loop on matrix rows

  // Now symmetrize the matrix:
  for (i = 0; i < BdkIsoOp::NOPS; ++i) { // loop on matrix rows
    for (int j = 0; j < i; ++j) {  // on matrix columns
      _evErr[i][j] = _evErr[j][i];
    }
  }
}
	      


const BdkTriBin & BdkBinning::bin(const BdkP2 & point) const {
  return _bins.at(binIndexContaining(point));
}

int BdkBinning::binIndexContaining(const BdkP2 & point) const {
  for (int b = 0; b < _bins.size(); ++b) {
    if (_bins[b].contains(point)) {
      return b;
    }
  }
  return -1;  // not found
}


int BdkBinning::binIndex(const BdkTriBin & bin) const {
  for (int b = 0; b < _bins.size(); ++b) {
    if (_bins[b] == bin) {  // operator== compares both center() and area() 
      return b;
    }
  }
  return -1;  // not found
}


 double BdkBinning::ev(int isoOp) const {
   assert(isoOp > -1 && isoOp < BdkIsoOp::NOPS);
   return _ev[isoOp];
 }



TCanvas * BdkBinning::Draw(Bool_t writeName,
			   Bool_t drawSymLines,
			   const char * drawDataOption) const {

  TCanvas * can = new TCanvas("Binning", "Binning", 800, 800);

  TH2F * h = new TH2F("Binning", "Binning", 30, 0, 
		      BdkDP::edge(), 30, 0, BdkDP::edge());
  
  h->SetXTitle("s_{12}");
  h->SetYTitle("s_{13}");
  h->Draw();

  if (TString(drawDataOption) == "SR" && 0 != _dataSR) {
    _dataSR->Draw("s13:s12", 0, "same");
  }
  else if (TString(drawDataOption) == "SB" && 0 != _dataSB) {
    _dataSB->Draw("s13:s12", 0, "same");
  }


  // Draw the Dalitz plot boundary:
  BdkDP::drawBoundary();
  
  if (drawSymLines) {
    // Draw the symmetry axes:
    BdkDP::drawE12Line();
    BdkDP::drawE13Line();
    BdkDP::drawE23Line();
  }

  // Draw all the bins:
  for (int b = 0; b < _bins.size(); ++b) {
    _bins[b].Draw(kBlack, 1, kFALSE, writeName);
  }

  return can;
}


// Tag a bin and its symmetry-related bins:
void BdkBinning::tagBin(int index, const char * nameChar) {
  TString name(nameChar);

  BdkTriBin & theBin = _bins.at(index);
  theBin.setName(name);
  
  BdkTriBin & binR = _bins[binIndex(theBin.R())];
  binR.setName(TString("R*") + name);
  
  BdkTriBin & binRR = _bins[binIndex(theBin.RR())];
  binRR.setName(TString("RR*") + name);
  
  BdkTriBin & binE = _bins[binIndex(theBin.E())];
  binE.setName(TString("E*") + name);
  
  BdkTriBin & binER = _bins[binIndex(theBin.ER())];
  binER.setName(TString("ER*") + name);
  
  BdkTriBin & binERR = _bins[binIndex(theBin.ERR())];
  binERR.setName(TString("ERR*") + name);
}


// Set all bin names to their nEventsSR():
void BdkBinning::setBinNamesNEventsSR() {
  for (int b = 0; b < _bins.size(); ++b) {
    BdkTriBin & bin = _bins[b];
    TString name;
    name += bin.nEventsSR();
    bin.setName(name);
  }
}


// Set all bin names to their serial number:
void BdkBinning::setBinNamesIndex() {
  for (int b = 0; b < _bins.size(); ++b) {
    BdkTriBin & bin = _bins[b];
    TString name;
    name += b;
    bin.setName(name);
  }
}


// Set all bin names to their nEventsSB():
void BdkBinning::setBinNamesNEventsSB() {
  for (int b = 0; b < _bins.size(); ++b) {
    BdkTriBin & bin = _bins[b];
    TString name;
    name += bin.nEventsSB();
    bin.setName(name);
  }
}


// Set all bin names to their inDalitzFrac():
void BdkBinning::setBinNamesInDalitzFrac() {
  for (int b = 0; b < _bins.size(); ++b) {
    BdkTriBin & bin = _bins[b];
    TString name;
    name += bin.inDalitzFrac();
    // clean the string of initial spaces:
    const char * nameStr = name;
    while(nameStr[0] == ' ') {
      ++nameStr;
    }
    bin.setName(nameStr);
  }
}

// Set all bin names to their waveFunc():
void BdkBinning::setBinNamesWaveFunc() {
  for (int b = 0; b < _bins.size(); ++b) {
    BdkTriBin & bin = _bins[b];
    TString name;
    name += bin.waveFunc();
    // clean the string of initial spaces:
    const char * nameStr = name;
    while(nameStr[0] == ' ') {
      ++nameStr;
    }
    bin.setName(nameStr);
  }
}

// Set all bin names to their waveFuncErr():
void BdkBinning::setBinNamesWaveFuncErr() {
  for (int b = 0; b < _bins.size(); ++b) {
    BdkTriBin & bin = _bins[b];
    TString name;
    name += bin.waveFuncErr();
    // clean the string of initial spaces:
    const char * nameStr = name;
    while(nameStr[0] == ' ') {
      ++nameStr;
    }
    bin.setName(nameStr);
  }
}


// Set all bin names to their efficiency():
void BdkBinning::setBinNamesEff() {
  for (int b = 0; b < _bins.size(); ++b) {
    BdkTriBin & bin = _bins[b];
    TString name;
    name += bin.efficiency();
    // clean the string of initial spaces:
    const char * nameStr = name;
    while(nameStr[0] == ' ') {
      ++nameStr;
    }
    bin.setName(nameStr);
  }
}


void BdkBinning::addBin(BdkTriBin & bin) {
  // Add the bin to the vector if it is at least partly inside the DP:
  BdkTriBin::InDalitzCode code = bin.inDalitzStat();
  if (BdkTriBin::NONE != code) {
    _bins.push_back(bin);
  }
}
