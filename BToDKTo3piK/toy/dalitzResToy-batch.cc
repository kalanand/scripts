void dalitzResToy_batch() {
  gROOT->ProcessLine(".x setup.cc");
  gROOT->ProcessLine(".L dalitzResToy.cc");
  setupPdf("dalitzRes-BadD.par");
}
  
