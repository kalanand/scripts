// calls fitChunks with progressive addition of bgds:

void fitProgressiveChunks(const char * directory = ".",
		 Bool_t replaceDeltaE = kFALSE,
		 const char * moreFixFile = "analysis/defaultMoreFixFile.cc",
		 const char * parInputFile = "analysis/defaultParInputFile.cc",
		 const char * floatFile = "analysis/defaultFloatFile.cc") {

  int nBins = 1.0 / WEIGHT_SIG;
  for (int b = 0; b < nBins; ++b) {
    TString fileName = directory;
    fileName += "/fitManyChunksig-";
    fileName += b;
    fileName += ".out";
 
    fitChunks(fileName.Data(), 
	      SIGBIT, 
	      b, 
	      replaceDeltaE, 
	      moreFixFile, 
	      parInputFile, 
	      floatFile);
  }
}
    

