#!/bin/tcsh -f

set resultDir = syst/mcStatSyst/

# include all par files into one big file
# no includes in sub par files allowed
cat ../BToDKTo3piK/params/all.par | grep include | awk '{print $2}' | xargs cat >! $resultDir/all.par

foreach t (`seq 0 9`)
    echo Submitting event type $t
    bsub -q xlong -J mcStatSyst-$t -o $resultDir/mcStatSyst-$t.log bbrroot -q -b setup.cc 'mcStatSyst.cc('$t',"'$resultDir'","'$resultDir/all.par'")'
end
