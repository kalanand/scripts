#!/bin/tcsh -f
# Run toy experiments with random CP parameters in rB/gamma/delta space
# Syntax: rndToy.csh N1 N2
#         where N1(2) is the number of the first(last) job
#

set nSamples = 50    # number of experiments per job
set dir = "/nfs/farm/babar/AWG17/BCK/Frank/toyMC/rndToy-2"  # output directory
#set dir = "/nfs/farm/babar/AWG17/BCK/Frank/toyMC/xyRndToy-2"  # output directory
#gROOT->ProcessLine(".x ../BToDKTo3piK/globals/setup.cc(BdkPdfDKDalitz::CART)");

set n1 = $1
set n2 = $2
foreach i (`seq $n1 $n2`)

# Create .cc file
cat <<EOF > $dir/fitMC-$i.cc
{
gROOT->ProcessLine(".x ../BToDKTo3piK/globals/setup.cc");
gROOT->ProcessLine(".L ../BToDKTo3piK/toy/rndToy.cc");
setupPdf();
BdkBatchMCStudy *toyMC = new BdkBatchMCStudy(pdfOnResDK,pdfOnResDK,RooArgSet(*m12,*m13,*Deltae,*qprime,*dprime),"e","mex",0,RooArgSet(*Hdtrkchge));
toyMC->useRandomCP();
toyMC->generateAndFitBatch("$dir/fitMC-$i.root",-1,$nSamples);
return 0;
}
EOF

# Submit the above root script to the batch queue
bsub -q kanga -C 0 -o $dir/fitMC-$i.log bbrroot -q -b $dir/fitMC-$i.cc

end
