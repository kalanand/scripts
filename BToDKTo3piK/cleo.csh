#!/usr/local/bin/tcsh -f
#
# This script runs GeneratorsQAApp once using the params in cleo.dec, then
# plots the result in paw. Overwrites last histogram file.
#

GeneratorsQAApp ../BToDKTo3piK/cleo.tcl
paw <<EOF 

hi/file 1  GfiMCTruth.hbook
opt zfl
opt nbox
zone 2 2
nt/pl 3.d0pppupmass*d0pppupmass%d0ppupmass*d0ppupmass
nt/pl 3.d0ppmupmass*d0ppmupmass%d0pppupmass*d0pppupmass
pict/print cleo.ps
EOF

gv cleo.ps &
