void signalToy_batch() {
  gROOT->ProcessLine(".x setup.cc");
  gROOT->ProcessLine(".L signalToy.cc");
  setupSimPdf();
}
  
