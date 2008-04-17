// Fit all available signal cunks:

void 
fitManyChunksSig2(int b1, int b2, const char * directory = "./",
		 const char * moreFixFile = "analysis/defaultMoreFixFile.cc",
		 const char * parInputFile = "analysis/defaultParInputFile.cc",
		 const char * floatFile = "analysis/defaultFloatFile.cc") {

  int nBins = 1.0 / WEIGHT_SIG ;
  int nbBins = 1.0 /WEIGHT_BCH ;
  int bb = 0;
  if( b2 > nBins ) {
	b2 = nBins;
  }

  for (int b = b1; b < b2; ++b) {
       TString fileName = directory;
       fileName += "fitManyChunksigbb-";
       fileName += b;
       fileName += "-";
       fileName += bb;
       fileName += ".out";
       int SIGBIT = SIG_G_BIT + SIG_B_BIT; 
       int AllBit = 2047;
       fitChunks2(fileName.Data(), 
	      255, 
	      b, bb, 
	      moreFixFile, 
	      parInputFile, 
	      floatFile);
       ++bb;
       if( bb >= int(1.0/WEIGHT_BCH) ) {
        	bb = 0;
       }
  }
}
    
