#! /usr/local/bin/tcsh -f
# Makes scripts to run fitShapedFlat in batch.

set xMStart = -0.022
set xMStep = 0.063
set yMStart = -0.019
set yMStep = 0.075


@ iXM = -1
while ($iXM < 2)
@ iXM = $iXM + 1

@ iYM = -1
while ($iYM < 2)
@ iYM = $iYM + 1

set fileBase = fitShapedFlat-$iXM-$iYM
set logFile =  $fileBase.log
set logFile2 =  $fileBase.log2
set datFile = $fileBase.dat
set file = $fileBase.csh

rm -f $file
rm -f $datFile
rm -f $logFile
rm -f $logFile2

echo "#! /usr/local/bin/tcsh -f" > $file
echo "date" >> $file
echo "bbrroot -l -b <<EOF >& $logFile2" >> $file
echo ".x ../BToDKTo3piK/globals/setup.cc " >> $file

set quotes = \"
echo "fitShapedFlat($quotes$datFile$quotes, $xMStart + $iXM * $xMStep, 0, 1, $yMStart + $iYM * $yMStep, 0, 1) " >> $file

echo ".q " >> $file
echo "EOF" >> $file
echo "date" >> $file

chmod +x $file

bsub -q long  $file -o $logFile

end
end
