#!/bin/tcsh -f

set parFile = "../BToDKTo3piK/params/syst/dstar_PPP_Group3.Bdk.par"
set title = modelToySyst
set resultDir = awg/toyMC/modelToySyst
set nExp = 20

foreach job (`seq 1 50`)
    echo Submitting job $job
    bsub -q xlong -J resAmpSyst_$job -o $resultDir/$title"_"$job.log bbrroot -q -b setup.cc 'systToyAnyParFile.cc("'$parFile'",'$nExp','$job',"'$title'","'$resultDir'")'
end
