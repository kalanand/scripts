#! /bin/tcsh
#
# sets up the scripts to run CLTreeMaker::makeCLTree and submits them.
#

set currentDir = `pwd`
set currentDir = `echo $currentDir:t`
if ("workdir" != $currentDir) then
  echo "This script should be run from workdir only. Exiting."
  exit
endif

echo -n "Enter number of jobs: "
set nJobs = $<

echo -n "Enter first job number, <cr> for 1: "
set firstJob = $<
if ("" == $firstJob) then
  set firstJob = 1
endif

echo -n "Enter number of points per job, <cr> for 1.e06: "
set nJobPoints = $<
if ("" == $nJobPoints) then
  set nJobPoints = 1000000
endif

echo -n "Enter bias factor: 1 for full bias, <cr> or 0 for none: "
set bias = $<
if ("" == $bias) then
  set bias = 0
endif

set dollar = \$

echo -n "Enter directory for csh and root files wrt. workdir,<cr> for ${dollar}MYWORK: "
set dir = $<
if ("" == $dir) then
  set dir = $MYWORK
endif


set dir = `echo $dir:r`
set fileBase = $dir/makeCLTree

echo -n "Enter queue, <CR> not to submit: "
set queue = $<

set quote = \"

@ job = 0
while ($job < $nJobs)
@ job = $job + 1
@ jobNum = $job + $firstJob - 1
set cshFile = $dir/makeCLTree-$jobNum.csh
rm -f $cshFile

echo "#! /bin/tcsh" > $cshFile
echo "bbrroot -l -b <<EOF" >> $cshFile
echo "CLTreeMaker maker;" >> $cshFile
echo "maker.makeCLTree($jobNum, $nJobPoints, $bias, $quote$fileBase$quote);" >> $cshFile
echo ".q" >> $cshFile
echo "EOF" >> $cshFile
chmod +x $cshFile

echo Wrote file $cshFile

if ($queue != "") then
  bsub -q $queue -o $MYWORK/makeCLTreeToy-$jobNum.log $cshFile
endif
end


