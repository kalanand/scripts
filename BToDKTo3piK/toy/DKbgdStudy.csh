#! /usr/local/bin/tcsh -f
# Makes scripts to run DKbgdStudy.cc in batch.

set quotes = \"

echo -n "Enter number of jobs: " 
set nJobs = $<

echo -n "Enter y for no bgd other than DK_bgd: "
set answer = $<
set noOtherBgd = 'kFALSE'
if ('y' == $answer) then
  set noOtherBgd = 'kTRUE'
endif

echo -n "Enter name of additional par file, 0 if none: "
set extraParFile = $<
if (0 != $extraParFile) then
  set extraParFile = $quotes$extraParFile$quotes  
endif

@ nJ = 0
while ($nJ < $nJobs)
@ nJ = $nJ + 1

set fileBase = DKbgdStudy-$nJ
set logFile =  $fileBase.log
set logFile2 =  $fileBase.log2
set outFile = $fileBase.dat
set file = $fileBase.csh

rm -f $file
rm -f $outFile
rm -f $logFile
rm -f $logFile2

echo "date" > $file
echo "bbrroot -l -b <<EOF >& $logFile2" >> $file
echo ".x ../BToDKTo3piK/globals/setup.cc " >> $file
echo "setRandomGenSeed($nJ)" >> $file

echo ".x ../BToDKTo3piK/toy/DKbgdStudy.cc($quotes$outFile$quotes, $noOtherBgd, kFALSE, 4, $extraParFile) " >> $file

echo ".q " >> $file
echo "EOF" >> $file
echo "date" >> $file

chmod +x $file

bsub -q medium  $file -o $logFile

end
