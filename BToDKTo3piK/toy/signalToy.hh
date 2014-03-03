//$Id: signalToy.hh,v 1.6 2006/03/20 19:05:54 fwinkl Exp $

BdkPdfDKDalitz *pdfBp;
BdkPdfDKDalitz *pdfBm;

// BdkBatchMCStudy doesn't know how to deal with pointers
RooCategory kcharge(*Hdtrkchge,"kcharge");

// Sim PDF
RooSimultaneous simPdf("simPdf","simPdf",kcharge);

// prototype dataset for Hdtrkchge
RooDataSet *proto;

