#! /bin/tcsh
#
# sets up the scripts to run CLTreeMaker::makeCLHist3 and submits them.
#

set currentDir = `pwd`
set currentDir = `echo $currentDir:t`
if ("workdir" != $currentDir) then
  echo "This script should be run from workdir only. Exiting."
  exit
endif

echo -n "Enter number of jobs for dividing rB, <cr> for 20: "
set nJobsIn = $<
if ("" == $nJobsIn) then
@ nJobs = 20
else
@ nJobs = $nJobsIn
endif

echo -n "Enter number of bins per variable per job, <cr> for 500: "
set nBinsGDIn = $<
if ("" == $nBinsGDIn) then
@ nBinsGD = 500
else
@ nBinsGD = $nBinsGDIn
endif

# nBins in rB is calculated so that after all jobs, it will be 
# the same asnBins in delta and gamma:
@ nBinsR = $nBinsGD / $nJobs

# make sure it looks like a double:
set nBinsR = $nBinsR.0

echo -n "Enter min rB for the entire set, <cr> for 0: "
set rBAllMin = $<
if ("" == $rBAllMin) then
  set rBAllMin = 0
endif

echo -n "Enter max rB for the entire set, <cr> for 1.0: "
set rBAllMax = $<
if ("" == $rBAllMax) then
  set rBAllMax = 1.0
endif


set dollar = \$

echo -n "Enter directory for csh and root files wrt. workdir,<cr> for ${dollar}MYWORK: "
set dir = $<
if ("" == $dir) then
  set dir = $MYWORK
endif


set dir = `echo $dir:r`
set fileBase = $dir/makeCLHist3

echo -n "Enter queue, <CR> not to submit (10M bins/job requires medium queue) : "
set queue = $<

set quote = \"

@ job = 0
while ($job < $nJobs)
@ job = $job + 1

# rBMin and rBMax are calculated to give the same density in rB as in the 
# angles, after all jobs are added together. 
set rBMin =  "($job - 1) * ($rBAllMax - $rBAllMin) * $nBinsR / $nBinsGD"
set rBMax =  "$job * ($rBAllMax - $rBAllMin) * $nBinsR / $nBinsGD"

set cshFile = $dir/makeCLHist3-$job.csh
rm -f $cshFile

echo "#! /bin/tcsh" > $cshFile
echo "bbrroot -l -b <<EOF" >> $cshFile
echo "CLTreeMaker maker;" >> $cshFile
echo "maker.makeCLHist3($nBinsR, $nBinsGD , $quote$fileBase-$job$quote, $rBMin, $rBMax, $rBAllMax);" >> $cshFile
echo ".q" >> $cshFile
echo "EOF" >> $cshFile
chmod +x $cshFile

echo Wrote file $cshFile

if ($queue != "") then
  bsub -q $queue -o $MYWORK/makeCLHist3-$job.log $cshFile
endif
end


