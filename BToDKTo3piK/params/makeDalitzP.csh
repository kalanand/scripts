#!/bin/tcsh -f
# Creates dalitzP.par from dalitzN.par
#
cat dalitzN.par | sed "s/dalitzHolderN/dalitzHolderP/" > dalitzP.par

