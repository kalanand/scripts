// Finds the #'s of events of each type in the MC, scaled to luminosity

void numEvtsMC(double lumiData = 112, 
	       double lumiOff = 99999, 
	       double lumiDK = 1148,
	       double lumiDpi = 99999,
	       double lumiB0 = 475,
	       double lumiBch = 468.4,
	       double lumiCc = 161.5,
	       double lumiUds = 166.8) {
  
  // Wake it up. For some reason this prevents crashes:
  data = read(chainDK);

  // DK good D:
  readCut = cutGoodD;
  data = read(chainDK);
  cout << "DK good D numEvts: ";
  if (0 != data) {
    cout << data->numEntries() 
	 << " normalized numEvts: " << data->numEntries() * lumiData / lumiDK;
  }
  cout << endl;

  
  // DK bad D:
  readCut = cutBadD;
  data = read(chainDK);
  cout << "DK bad D numEvts: ";
  if (0 != data) {
    cout << data->numEntries() 
	 << " normalized numEvts: " << data->numEntries() * lumiData / lumiDK;
  }
  cout << endl;

  
  
  // Dpi good D:
  readCut = cutGoodD;
  data = read(chainDpi);
  cout << "Dpi good D numEvts: ";
  if (0 != data) {
    cout << data->numEntries() 
	 << " normalized numEvts: " << data->numEntries() * lumiData / lumiDpi;
  }
  cout << endl;

  
  // Dpi bad D:
  readCut = cutBadD;
  data = read(chainDK);
  cout << "Dpi bad D numEvts: ";
  if (0 != data) {
    cout << data->numEntries() 
	 << " normalized numEvts: " << data->numEntries() * lumiData / lumiDpi;
  }
  cout << endl;

  
  
  // charmless:
  readCut = "";
  data = read(chainCharmless);
  cout << "charmless numEvts: ";
  if (0 != data) {
    cout << data->numEntries() 
	 << " normalized numEvts: " 
	 << data->numEntries() * lumiData / ((lumiB0 + lumiBch)/2;)
  }
  cout << endl;

  
  
  // Dpipi:
  readCut = cutGoodD;
  data = read(chainDpipi);
  cout << "Dpipi numEvts: ";
  if (0 != data) {
    cout << data->numEntries() 
	 << " normalized numEvts: " 
	 << data->numEntries() * lumiData / ((lumiB0 + lumiBch)/2;)
  }
  cout << endl;

  
  
  // B+ bad D:
  readCut = cutBadD;
  data = read(chainBchComb);
  cout << "B+ bad D numEvts: ";
  if (0 != data) {
    cout << data->numEntries() 
	 << " normalized numEvts: " << data->numEntries() * lumiData / lumiBch;
  }
  cout << endl;

  
  // B+ good D:
  readCut = cutGoodD;
  data = read(chainBchComb);
  cout << "B+ good D numEvts: ";
  if (0 != data) {
    cout << data->numEntries() 
	 << " normalized numEvts: " << data->numEntries() * lumiData / lumiBch;
  }
  cout << endl;

  
  
  // B0 bad D:
  readCut = cutBadD;
  data = read(chainB0Comb);
  cout << "B0 bad D numEvts: ";
  if (0 != data) {
    cout << data->numEntries() 
	 << " normalized numEvts: " << data->numEntries() * lumiData / lumiB0;
  }
  cout << endl;

  
  // B0 good D:  
  readCut = cutGoodD;
  data = read(chainB0Comb);
  cout << "B0 good D numEvts: ";
  if (0 != data) {
    cout << << data->numEntries() 
	 << " normalized numEvts: " << data->numEntries() * lumiData / lumiB0;
  }
  cout << endl;

  
  
  // uds:
  readCut = cutBadD;
  data = read(chainUds);
  cout << "uds numEvts: ";
  if (0 != data) {
    cout << data->numEntries() 
	 << " normalized numEvts: " << data->numEntries() * lumiData / lumiUds;
  }
  cout << endl;

  
  
  // cc bad D:
  readCut = cutBadD;
  data = read(chainCc);
  cout << "cc bad D numEvts: ";
  if (0 != data) {
    cout << data->numEntries() 
	 << " normalized numEvts: " << data->numEntries() * lumiData / lumiCc;
  }
  cout << endl;

  
  // cc good D:
  readCut = cutGoodD;
  data = read(chainCc);
  cout << "cc good D numEvts: ";
  if (0 != data) {
    cout << data->numEntries() 
	 << " normalized numEvts: " << data->numEntries() * lumiData / lumiCc;
  }
  cout << endl;

  
  
  // off resonance:
  readCut = "";
  data = read(chainOffData);
  cout << "off resonance numEvts: ";
  if (0 != data) {
    cout << data->numEntries() 
	 << " normalized numEvts: " << data->numEntries() * lumiData / lumiOff;
  }
  cout << endl;

}
