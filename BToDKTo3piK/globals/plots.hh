// Declares the global plots

#ifndef PLOTS_HH
#define PLOTS_HH


// frames:
RooPlot * mesFrame = 0;
RooPlot * DeltaeFrame = 0;
RooPlot * d0massFrame = 0;
RooPlot * nnFrame = 0;
RooPlot * bknnFrame = 0;
RooPlot * R2Frame = 0;
RooPlot * qprimeFrame = 0;
RooPlot * dprimeFrame = 0;
RooPlot * m12Frame = 0;
RooPlot * m13Frame = 0;

// Default file name:
TString plotFileName = "plots.eps";

// same plots with likelihood cuts:
RooPlot * mesFrameLike = 0;
RooPlot * DeltaeFrameLike = 0;
RooPlot * d0massFrameLike = 0;
RooPlot * nnFrameLike = 0;
RooPlot * bknnFrameLike = 0;

// Default file name:
TString plotFileNameLike = "plotsLike.eps";

#endif


