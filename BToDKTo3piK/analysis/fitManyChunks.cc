// call fitChunks many times

void fitManyChunks(const char * resultFileBaseName, 
		   int todo = 0, 
		   int startBin = 0,
		   int endBin = 0,
		   Replacer replacer = Replacer(), 
		   const char * moreFixFile = "analysis/defaultMoreFixFile.cc",
		   const char * parInputFile = "analysis/defaultParInputFile.cc",
		   const char * floatFile = "analysis/defaultFloatFile.cc") {
  
  doPlot = kFALSE;

  for (int bin = startBin; bin < endBin; ++bin) {
    TString fileName = TString(resultFileBaseName);
    fileName += "-";
    fileName += bin;
    
    fitChunks(fileName, todo, bin, replacer,
	      moreFixFile,
	      parInputFile,
	      floatFile);
  }
}
