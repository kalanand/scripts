// calls fitChunks with progressive addition of bgds:

// this function is just so as not to repeat code:
void 
runFitProgChunks(int todo, 
		 int bin, 
		 const char * directory = "./",
		 Replacer replacer = Replacer(), 
		 const char * moreFixFile = "analysis/defaultMoreFixFile.cc",
		 const char * parInputFile = "analysis/defaultParInputFile.cc",
		 const char * floatFile = "analysis/defaultFloatFile.cc") {
  
    TString baseFile = directory;
    baseFile += todo;
    baseFile += "-";
    baseFile += bin;
    outFile = baseFile + ".out";
    plotFileName = baseFile + ".eps";
    
    fitChunks(outFile.Data(), 
	      todo, 
	      bin, 
	      replacer, 
	      moreFixFile, 
	      parInputFile, 
	      floatFile);
};

//------------------------------------------------------------------
void 
fitProgChunks(const char * directory = "./",
	      int nBins = 1, // 0 for max # of bins
	      Replacer replacer = Replacer(), 
	      const char * moreFixFile = "analysis/defaultMoreFixFile.cc",
	      const char * parInputFile = "analysis/defaultParInputFile.cc",
	      const char * floatFile = "analysis/defaultFloatFile.cc") {
  
  if (0 == nBins) {
    nBins = 1.0 / WEIGHT_SIG;
  }

  for (int b = 0; b < nBins; ++b) {
    runFitProgChunks(SIG_G_BIT + SIG_B_BIT, 
		     b, directory, replacer, 
		     moreFixFile,
		     parInputFile,
		     floatFile);
    
    runFitProgChunks(SIG_G_BIT + BB_B_BIT, 
		     b, directory, replacer, 
		     moreFixFile,
		     parInputFile,
		     floatFile);
    
    runFitProgChunks(SIG_G_BIT + SIG_B_BIT + QQ_B_BIT, 
		     b, directory, replacer, 
		     moreFixFile,
		     parInputFile,
		     floatFile);
    
    runFitProgChunks(SIG_G_BIT + SIG_B_BIT + BB_B_BIT,
		     b, directory, replacer, 
		     moreFixFile,
		     parInputFile,
		     floatFile);


    runFitProgChunks(SIG_G_BIT + SIG_B_BIT + BB_B_BIT + QQ_B_BIT,
		     b, directory, replacer, 
		     moreFixFile,
		     parInputFile,
		     floatFile);


    runFitProgChunks(SIG_G_BIT + SIG_B_BIT + BB_B_BIT + QQ_B_BIT + QQ_G_BIT,
		     b, directory, replacer, 
		     moreFixFile,
		     parInputFile,
		     floatFile);


    runFitProgChunks(SIG_G_BIT + SIG_B_BIT + BB_B_BIT + QQ_B_BIT + QQ_G_BIT + DPI_G_BIT,
		     b, directory, replacer, 
		     moreFixFile,
		     parInputFile,
		     floatFile);



  }
}
