// Declares the PDFs and their parameters:

#ifndef PDFS_HH
#define PDFS_HH

#include "BToDKTo3piK/BdkPdfOnRes.hh"
#include "BToDKTo3piK/BdkPdfHistPdfHolder.hh"
#include "BToDKTo3piK/BdkDalitzPdfHolder1.hh"
#include "BToDKTo3piK/BdkPdf2DpolyDalitz.hh"
#include "BToDKTo3piK/BdkDalitzCfg.hh"
#include "BToDKTo3piK/BdkPdfDKDalitz.hh"
#include "BToDKTo3piK/BdkPdfDDalitzInc.hh"
#include "BToDKTo3piK/BdkPdfAbsBase.hh"

BdkPdfOnRes pdfOnResDK;  // for the DK sample

// holders for the NN histogram PDFs:
BdkPdfHistPdfHolder nnHolder;
BdkPdfHistPdfHolder bknnHolder;

// holders for Dalitz shapes for positive and negative B candidates:
BdkDalitzPdfHolder1 dalitzHolderP;
BdkDalitzPdfHolder1 dalitzHolderN;

// Efficiency functions
BdkPdf2DpolyDalitz eff;            // signal
BdkPdf2DpolyDalitz effOther;       // all other types
BdkPdf2DpolyDalitz effDKBadD;      // DKBadD efficiency (for systematics)

// Global Dalitz plot configuration (like m23 mass veto window)
BdkDalitzCfg *dalitzCfg;

// A D Pdf, and a signal and a background PDF to play with:
BdkPdfDDalitz dDalitz;
BdkPdfDKDalitz dkDalitz;
BdkPdfDDalitzInc bgdDalitz;

BdkPdfAbsBase* dalitzRes12;
BdkPdfAbsBase* dalitzRes13;

// The following pdfs are setup in setupAuxPdfs
BdkPdfDKDalitz* ksppPdf;
BdkDDalitzKsPiPiAmp* ksppAmp;

BdkPdfDKDalitz* kskkPdf;
BdkDDalitzKsKKAmp* kskkAmp;

BdkPdfDKDalitz* kskpPdf;
BdkPdfDDalitz* kskpDPdf;
BdkPdfDDalitz* kskpDbarPdf;
BdkDDalitzKsKPiAmp* kskpAmp;
BdkDDalitzKsKPiAmp* kskpBarAmp;

BdkPdfDKDalitz* kkpPdf;
BdkDDalitzKKPiAmp* kkpAmp;

#endif

