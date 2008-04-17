#!/bin/tcsh -f

set nEvents = 50000

bsub -q xlong bbrroot -q -b setup.cc setupAuxPdfs.cc 'plotDalitzAmp.cc(dalitzHolderN.sigGoodD0Type(),'$nEvents')'
bsub -q xlong bbrroot -q -b setup.cc setupAuxPdfs.cc 'plotDalitzAmp.cc(*ksppPdf,'$nEvents')'
bsub -q xlong bbrroot -q -b setup.cc setupAuxPdfs.cc 'plotDalitzAmp.cc(*kkpPdf,'$nEvents')'
bsub -q xlong bbrroot -q -b setup.cc setupAuxPdfs.cc 'plotDalitzAmp.cc(*kskkPdf,'$nEvents',"m12:m23")'
bsub -q xlong bbrroot -q -b setup.cc setupAuxPdfs.cc 'plotDalitzAmp.cc(*kskpDPdf,'$nEvents',"m23:m12","kskpDPdf")'
bsub -q xlong bbrroot -q -b setup.cc setupAuxPdfs.cc 'plotDalitzAmp.cc(*kskpDbarPdf,'$nEvents',"m23:m12","kskpDbarPdf")'
