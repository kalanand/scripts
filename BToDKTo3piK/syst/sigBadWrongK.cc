// $Id: sigBadWrongK.cc,v 1.1 2006/09/25 23:31:13 fwinkl Exp $
// Study the effect of wrong sign kaons in DKBadD events
//
// Replace Dalitz variables in DKBadD events with a bad K+
// with toy events from B- signal and vice versa.
// Get a new DKBadD Dalitz histogram.
// Requires getDalitzParams.cc being loaded.

void sigBadWrongK()
{
  // read DKBadD events with a good D (hence bad K)
  readCut = cutDKBadD+cutGoodD;
  RooDataSet* dataBadK = read(sigTree);

  // Replace m12/m13 for events with a K+ with toy events 
  // generated from the B- signal PDF and vice versa.
  // Notice that replace() expects the B+ PDF as the first argument
  // We give the B- PDF as first argument to do the flipping
  RooDataSet* dataGoodK = replace(dataBadK,RooArgSet(*m12,*m13),
				  dalitzHolderN.sigGoodD0Type(),
				  dalitzHolderP.sigGoodD0Type());

  // now read DKBadD events with a bad D (hence good K)
  readCut = cutDKBadD+cutBadD;
  data = read(sigTree);

  // merge them with the toy events and get new Dalitz shape
  data->append(*dataGoodK);
  getDalitzParamsDKBadD(&data->tree(),"");
}
