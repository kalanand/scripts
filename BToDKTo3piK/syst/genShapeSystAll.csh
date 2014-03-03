#
# run genShapeSyst for all par files in params/DK/
#

set files = "deltaE.par"
foreach file ($files)
  set f = `echo $file:t`
  rm -f genShapeSyst-$f.log
  echo running on $f > genShapeSyst-$f.log
  genShapeSyst <<EOF >>& genShapeSyst-$f.log
params/DK/
syst/
$f
EOF
end
