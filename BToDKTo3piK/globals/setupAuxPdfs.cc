// $Id: setupAuxPdfs.cc,v 1.3 2008/01/02 20:46:06 fwinkl Exp $
// Setup auxiliary pdfs for additional studies

#ifndef SETUPAUXPDFS_CC
#define SETUPAUXPDFS_CC

void setupAuxPdfs(BdkPdfDKDalitz::COORD coord = BdkPdfDKDalitz::POLAR)
{
  // Remove veto window (this is global!)
  dalitzCfg->removeM23Veto();
  
  //
  // pdfs and amplitdues for D0 -> Ks(1) pi+(2) pi-(3)
  //
  ksppAmp = new BdkDDalitzKsPiPiAmp("ksppAmp","Ks pi+ pi-",0);
  ksppPdf = new BdkPdfDKDalitz("ksppPdf","Ks pi+ pi-",*m12,*m13,
                               BdkDalitzBase::D0,coord, 0, ksppAmp,
                               BdkDalitz::PPK0);

  // Remove efficiency function
  // not really needed because efficiency() is overwritten anyways to return 1.0
  ksppPdf->setEfficiencyFunc(0);
  
  //  ksppPdf->setVerbose("dgciv+");
  ksppAmp->setPdf(ksppPdf->pdfType());
  ksppPdf->parameters().readFromFile("../BToDKTo3piK/params/ksppPdf.par");
  
  // Flip phases of KsPi- spin 1 resonances (due to different def. of spin factor)
  // You need to recalculate x0 after this !!!
  for (int i=0; i<ksppAmp->nComps(); i++) {   
    if (ksppAmp->trackinfo(i)==BdkAbsDDalitzAmp::PIM_PI0 && ksppAmp->spinRes(i)==1)
      ksppAmp->phaseRes(i)->setVal(ksppAmp->phaseRes(i)->getVal()-180);    
  }
  // ksppPdf->recalcX0();

  //
  // pdfs and amplitudes for D0 -> pi0(1) K+(2) K-(3)
  //
  kkpAmp = new BdkDDalitzKKPiAmp("kkpAmp","pi0 K+ K-",0);
  kkpPdf = new BdkPdfDKDalitz("kkpPdf","pi0 K+ K-",*m12,*m13,
                              BdkDalitzBase::D0,coord, 0, kkpAmp,
                              BdkDalitz::KKP0);

  kkpPdf->setEfficiencyFunc(0);
  kkpAmp->setPdf(kkpPdf->pdfType());
  kkpPdf->parameters().readFromFile("../BToDKTo3piK/params/kkpPdf.par");

  //
  // pdfs and amplitdues for D0 -> Ks(1) K+(2) K-(3)
  //
  kskkAmp = new BdkDDalitzKsKKAmp("kskkAmp","Ks K+ K-",0);
  kskkPdf = new BdkPdfDKDalitz("kskkPdf","Ks K+ K-",*m12,*m13,
                               BdkDalitzBase::D0,coord, 0, kskkAmp,
                               BdkDalitz::KKK0);

  kskkPdf->setEfficiencyFunc(0);
  kskkAmp->setPdf(kskkPdf->pdfType());
  kskkPdf->parameters().readFromFile("../BToDKTo3piK/params/kskkPdf.par");  
  
  //
  // pdfs and amplitdues for D0 -> Ks(1) pi+(2) K-(3)
  //
  // This PDF uses two separate D/Dbar amplitudes
  // To normalize everything properly do:
  //    kskpAmp->calNorm();
  //    kskpBarAmp->calNorm();
  //    ((BdkDKNonCDalitz*)kskpPdf->getPdf())->calDDbarNorm();
  
  kskpAmp = new BdkDDalitzKsKPiAmp("kskpAmp","D0 -> Ks pi+ K-",0,
                                   BdkDDalitzKsKPiAmp::allComponents(),+1);
  kskpBarAmp = new BdkDDalitzKsKPiAmp("kskpBarAmp","D0bar -> Ks pi+ K-",0,
                                      BdkDDalitzKsKPiAmp::allComponents(),+1);

  
  kskpPdf = new BdkPdfDKDalitz("kskpPdf","Ks pi+ K-",*m12,*m13,
                               BdkDalitzBase::D0,coord, 0,
                               kskpAmp, BdkDalitz::KPK0, kskpBarAmp);

  kskpDPdf = new BdkPdfDDalitz("kskpDPdf","D0 -> Ks pi+ K-",*m12,*m13,
                               BdkDalitzBase::D0, kskpAmp, BdkDDalitzKsKPiAmp::allComponents(),
                               1, BdkDalitz::KPK0);

  // We don't want this PDF to flip m12/m13 itself. Therefore needs to be D0 flavor!
  kskpDbarPdf = new BdkPdfDDalitz("kskpDbarPdf","D0bar -> Ks pi+ K-",*m12,*m13,
                                  BdkDalitzBase::D0, kskpBarAmp, BdkDDalitzKsKPiAmp::allComponents(),
                                  1, BdkDalitz::KPK0);

  kskpPdf->setEfficiencyFunc(0);
  kskpAmp->setPdf(kskpPdf->pdfType());
  kskpBarAmp->setPdf(kskpPdf->pdfType());
  kskpPdf->parameters().readFromFile("../BToDKTo3piK/params/kskpPdf.par");
}

#endif
