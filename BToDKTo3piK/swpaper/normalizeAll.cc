// $Id: normalizeAll.cc,v 1.2 2008/01/03 18:11:23 fwinkl Exp $
// Normalize all auxiliary PDFs
// Use the normalizeAll.sh script to merge the new normalization constants
// with the original par files in params/

void normalizeAll(int events = (int)1e6)
{
  normalize(ksppPdf,events);  
  normalize(kkpPdf,events);
  normalize(kskkPdf,events);
  normalizeKsKPi(events);
}

void normalize(BdkPdfDKDalitz* pdf, int events)
{
  pdf->dalitzAmp()->calNorm(events);
  // Update x0/y0
  pdf->x0()->setVal(pdf->dalitzAmp()->x0());
  pdf->y0()->setVal(pdf->dalitzAmp()->y0());
  pdf->parameters().writeToFile(pdf->GetName());
}

void normalizeKsKPi(int events)
{
  // Special for KsKPi since two separate D Dalitz amps
  kskpAmp->calNorm(events);
  kskpBarAmp->calNorm(events);
  ((BdkDKNonCDalitz*)kskpPdf->getPdf())->calDDbarNorm(events);
  kskpPdf->x0()->setVal(kskpAmp->x0());
  kskpPdf->y0()->setVal(kskpAmp->y0());
  kskpPdf->parameters().writeToFile(kskpPdf->GetName());
}
