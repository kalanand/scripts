// $Id: countEvtTypes.cc,v 1.5 2006/04/25 17:51:47 fwinkl Exp $
// Simple script to count the different event types

void countEvtTypes(TTree* sigTree = sigTree)
{
  TCut baseCut = cutSigReg;

  // This list has the following format:
  // tree1, cut1, cut2, cut3, tree2, cut4, cut5, tree3, cut6, ... cutN
  
  TList list;

  list.Add(bpTree);
  list.Add(&cutBBBadD);
  list.Add(&cutBBGoodD);
  list.Add(&cutDPiX);
  list.Add(&cutDKX);

  list.Add(b0Tree);
  list.Add(&cutBBBadD);
  list.Add(&cutBBGoodD);
  list.Add(&cutDPiX);
  list.Add(&cutDKX);
  
  list.Add(bbTree);
  list.Add(&cutBBBadD);
  list.Add(&cutBBGoodD);
  list.Add(&cutDPiX);
  list.Add(&cutDKX);
  
  list.Add(dpiTree);
  list.Add(&cutDPiBadD);
  list.Add(&cutDPiGoodD);
  
  list.Add(ccTree);
  list.Add(&cutqqBadD);
  list.Add(&cutqqGoodD);

  list.Add(udsTree);
  list.Add(&cutqqBadD);
  list.Add(&cutqqGoodD);

  list.Add(qqTree);
  list.Add(&cutqqBadD);
  list.Add(&cutqqGoodD);
  
  gROOT->cd();

  cout <<"Expected number of events based on "<<ON_PEAK_LUMI<<" fb-1 of on-peak data."
       <<endl<<endl;
  
  printf("%-10s%-10s%6s%10s%10s\n","Type","MC","events","expected","error");
  int i = 0;
  while (i<list.GetSize()) {
    // Apply base cut to tree
    TTree *tree = ((TTree*)list.At(i))->CopyTree(baseCut);
    tree->SetTitle(list.At(i)->GetTitle());
    // Get the sample weight
    Double_t weight = getTreeWeight((TTree*)list.At(i));

    // Loop over cuts
    while (++i<list.GetSize() && TString(list.At(i)->ClassName())=="TCut") {
      TCut* cut = (TCut*)list.At(i);
      Int_t events = (Int_t)tree->Draw("mes",*cut,"goff");
      Double_t expected = events*weight;
      Double_t error = sqrt(events)*weight;
        
      printf("%-10s%-10s%6.f%10.1f%10.1f\n",
	     cut->GetName(),tree->GetTitle(),events,expected,error);
    }
    delete tree;
  }
}

