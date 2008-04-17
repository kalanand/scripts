// $Id: drawBinning.cc,v 1.1 2006/02/18 23:31:18 fwinkl Exp $
// Draw the bin boundaries on a 2D histogram

void drawBinning(TH2 *h) {

  if (!h) return;
  if (h->GetDimension()!=2) return;

  TLine line;
  TAxis *ax = h->GetXaxis();
  TAxis *ay = h->GetYaxis();
  if (!ax || !ay) return;

  for (Int_t i=ax->GetFirst()+1; i<=ax->GetNbins(); i++) {
    line.DrawLine(ax->GetBinLowEdge(i),ay->GetBinLowEdge(ay->GetFirst()),
                  ax->GetBinLowEdge(i),ay->GetBinUpEdge(ay->GetLast()));
  }

  for (Int_t i=ay->GetFirst()+1; i<=ay->GetNbins(); i++) {
    line.DrawLine(ax->GetBinLowEdge(ax->GetFirst()),ay->GetBinLowEdge(i),
                  ax->GetBinUpEdge(ax->GetLast()),ay->GetBinLowEdge(i));
  }
}
