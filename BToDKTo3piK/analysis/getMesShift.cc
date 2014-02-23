void getMesShift(RooAbsData *dat) {
	BdkPdfCBArgus * mesPdf = new BdkPdfCBArgus("mesFit", "mes fit", *mes);
        mesPdf->parameters()->readFromFile("mesData.txt");
	RooFitResult * Result =  mesPdf->getPdf()->fitTo(*dat);
	if( 0 != Result ) {
	    Result->Print("V");
	}
}
 

