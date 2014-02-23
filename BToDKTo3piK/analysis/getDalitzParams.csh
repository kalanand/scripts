#! /usr/local/bin/tcsh -f
# $Id: getDalitzParams.csh,v 1.2 2006/04/25 17:51:48 fwinkl Exp $
# Runs the getDalitzParams root script
#

rm -f getDalitzParams.log
bbrroot -b -l <<EOF >& getDalitzParams.log
.x ../BToDKTo3piK/globals/setup.cc
.L ../BToDKTo3piK/analysis/getDalitzParams.cc
getDalitzParamsAll(ALL_TYPES,true,true)
.q
EOF
