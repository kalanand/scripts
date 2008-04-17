#! /usr/local/bin/tcsh -f
# because we're having trouble running getParams for all types without 
# crashing, run it in many jobs:
#

#set vars = "d0mass Deltae mes qprime"
set vars = "Deltae qprime dprime"

foreach var ($vars)
@ n = 0
#while ($n < 11)
#@ n = $n + 1
rm -f getParams-$var.log
bbrroot -b -l <<EOF >& getParams-$var.log
.x ../BToDKTo3piK/globals/setup.cc
// doFit = kFALSE
getParamsAll($var, ALL_TYPES, true)
.q
EOF
#end
end
