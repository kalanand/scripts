// Fit all available signal cunks:

void 
fitManyChunksSig(const char * directory = "./",
		 Bool_t replaceDeltaE = kFALSE,
		 const char * moreFixFile = "analysis/defaultMoreFixFile.cc",
		 const char * parInputFile = "analysis/defaultParInputFile.cc",
		 const char * floatFile = "analysis/defaultFloatFile.cc") {

  int nBins = 1.0 / WEIGHT_SIG ;
  for (int b = 50; b < 90 ; ++b) {
    TString fileName = directory;
    fileName += "fitManySigBB-";
    fileName += b;
    fileName += ".out";
    int SIGBIT = SIG_G_BIT + SIG_B_BIT; 
    int AllBit = 2047;
    int AllBB = 255;
    fitChunks(fileName.Data(), 
	      AllBB, 
	      b, 
	      replaceDeltaE, 
	      moreFixFile, 
	      parInputFile, 
	      floatFile);
  }
}
    
