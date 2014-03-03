// $Id: mergeFiles.cc,v 1.6 2006/05/29 22:57:14 fwinkl Exp $
// merge R16 and R18 MC root files

void mergeR16R18() {

  //  TString outDir = "/nfs/farm/babar/AWG17/BCK/Frank/ntuples/R16R18-MC/";
  TString outDir = "/nfs/farm/babar/AWGbreco02/BDK/GoodSmall/SPmerge_F/";
  merge2Files(sigFlatFile16, sigFlatFile18, outDir);
  merge2Files(sigFile16, sigFile18, outDir);
  merge2Files(bbFile16, bbFile18, outDir);
  merge2Files(qqFile16, qqFile18, outDir);
  merge2Files(dpiFile16, dpiFile18, outDir);

  // for DPI
  //  merge2Files(dpiBbFile16,dpiBbFile18,"GoodSmall/dpi_F/SPmerge/");
  //  merge2Files(dpiQqFile16,dpiQqFile18,"GoodSmall/dpi_F/SPmerge/");
  //  merge2Files(dpiDpiFile16,dpiDpiFile18,"GoodSmall/dpi_F/SPmerge/");

}

void merge2Files(TFile* f1, TFile* f2, TString outDir)
{
  if (f1==0 || f2==0) {
    cout << "ERROR: One of the two files does not exist."<<endl;
    return;
  }

  TChain c("h1");
  c.Add(f1->GetName());
  c.Add(f2->GetName());

  // Need to close the original files otherwise the merge doesn't work
  f1->Close();
  f2->Close();
  Long64_t N = c.GetEntries();
  cout << "Merging "<<N<<" entries from"<<endl;
  c.ls();
  TString filename = outDir+TString(f1->GetName()).Tokenize("/")->Last()->GetName();
  cout << "into "<<filename<<endl;
  c.Merge(filename);

  // Check the merged file
  TFile f(filename);
  TTree* tree = (TTree*)f.Get("h1");
  Long64_t N2 = tree->GetEntries();
  if (N!=N2) cout << "ERROR: Number of events do not match"<<N<<"!="<<N2<<endl;
  f.Close();
}
