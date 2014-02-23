#!/bin/tcsh -f 
# Shell script to run fitMC.cc

set file = "awg/fitMC/fitMC-1"

# Make a copy of script
cp ../BToDKTo3piK/analysis/fitMC.cc $file.cc
rm -f $file.log
rm -f $file.root


# DPIX_BIT + BB_B_BIT
# DPIX_BIT + BB_B_BIT + DKX_BIT
# -QQ_G_BIT-QQ_B_BIT
# ALL_TYPES
# SIG_G_BIT+SIG_B_BIT+DPI_B_BIT+DKX_BIT+BB_B_BIT+DPIX_BIT+DPI_G_BIT

# new RooArgSet(*m12,*m13),DPI_G_BIT,
# compFitMC(ALL_TYPES-QQ_G_BIT-BB_G_BIT,"$file",-1,new RooArgSet(*m12,*m13),DPI_G_BIT,5.0);
# compFitMC(ALL_TYPES-QQ_G_BIT-QQ_B_BIT,"$file",-1,0,0,5.0);

cat <<EOF > $file-run.cc
{
gROOT->ProcessLine(".L ../BToDKTo3piK/utils/TypeBits.hh");
gROOT->ProcessLine(".x setupVars.cc");
gSystem->CompileMacro("$file.cc","fg");
compFitMC(ALL_TYPES-QQ_G_BIT-BB_G_BIT,"$file",-1,new RooArgSet(*m12,*m13),DPI_G_BIT,1.0);
}
EOF

bsub -q kanga -C 0 -R cob -o $file.log -J $file bbrroot -q -b $file-run.cc


# (polar toys)
# 1   all types, replace DpiGoodD Dalitz, remove BBgood, qqGood
# 2   all types, replace DpiGoodD Dalitz, remove BBgood, qqGood, scale x5
# 3   no qqbar events, scale x5


# in old4/
#
# (cartesian toys)
# 1   all types, replace DpiGoodD Dalitz, remove BBgood, qqGood
# 2   all types, replace DpiGoodD Dalitz, remove BBgood, qqGood, scale x5
# 3   no qqbar events, scale x5



# in old3/
#
# 1   all types
# 2   no qqbar events, scale x3
# 3   no qqbar events, replace DpiGoodD Dalitz with DDalitzInc PDF, scale x3
# 4   no qqbar events, replace DpiGoodD Dalitz, scale x3
# 5   no qqbar events, scale x5
# 6   all types, replace DpiGoodD Dalitz
# 7   no qqbar events, replace DpiGoodD Dalitz with DDalitzInc PDF, scale x5
# 8   no qqbar events, replace DpiGoodD Dalitz, scale x5
# 9   all types, replace DpiGoodD Dalitz, remove BBgood, qqGood, scale x5
# 10  all types, replace DpiGoodD Dalitz, remove BBgood, qqGood
# 11  all types, no Dalitz
# 12  all types, x/y BLINDED




# in old2/
#
# 1   all types
# 2   all types, no signal
# 3   all types, replace Dalitz for DPiX
# 4   all types, no signal, replace Dalitz for DPiX
# 5   all types, no signal, replace qprime for DPiX
# 6   no qqbar events
# 7   signal only
# 8   DPiX, DKX
# 9   DPiX, DKX, BBBadD
# 10  DPiX, DKX, BBBadD, DpiBadD
# 11  DPiX, DKX, qqBadD
# 12  DPiX, qqBadD
# 13  qqBadD, BBBadD

# 14  no qqbar events, scale x 3
# 15  no qqbar, signal, no good D background, scale x 3
# 16  = 15, but no DPiX
# 17  no qqbar, signal, DpiGoodD, badD background, scale x 3
# 18  = 17, replace DpiGoodD Dalitz
# 19  = 17, replace DpiGoodD Dalitz with DDalitzInc PDF





# in old/
#
# 11 signal eff.
# 12 = 11 with background eff.
# 13 = 12 only replace Dalitz for DPIX (bkg eff)
# 14 = 12 only replace qprime (bkg eff)
# 15 all types, no signal, bkg eff. replace Dalitz for DPIX
# 16 all types, no signal, bkg eff.
# 17 all types, no signal, sig eff. replace Dalitz for DPIX
# 18 all types, no signal, sig eff.
