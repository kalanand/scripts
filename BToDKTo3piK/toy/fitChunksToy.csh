#! /usr/local/bin/tcsh -f
#
# This script runs or submits jobs of fitToy.cc in a loop
#

echo -n "Enter directory for results (.out) files: "
set dir = $<

if ("" == "$dir") then
  set dir = "."
endif

if (!(-d $dir)) then
  echo "Directory $dir does not exist. Making it."
  mkdir $dir
endif

echo -n "Enter first job number: " 
set firstJob =  $<

echo -n "Enter last job number: " 
set lastJob =  $<

echo -n "Enter the start counter: "
set stc = $<

echo -n "Enter the end counter: "
set etc = $<

echo -n "Which queue do you want to submit the job: "
set queue = $<

if ("" == $firstJob || "" == $lastJob) then
  echo empty argument given for either firstJob, lastJob. Exiting.
  exit
endif

if("" == $stc || "" == $etc ) then
  echo empty argument given for either first counter and last counter. Exiting.
  exit
endif 

echo -n "Enter parInputFile (path from workdir), <return> for default: "
set parInputFile = $<

if ("" == $parInputFile) then
  set parInputFile = analysis/defaultParInputFile.cc
endif
if (!(-e $parInputFile)) then
  echo "$parInputFile does not exist. Exiting."
  exit
endif

echo -n "Enter floatFile (path from workdir), <return> for default: "
set floatFile = $<

if ("" == $floatFile)then
  set floatFile = analysis/defaultFloatFile.cc
endif
if (!(-e $floatFile)) then
  echo "$floatFile does not exist. Exiting."
  exit
endif

echo -n "Enter fixFile (path from workdir), <return> for default: "
set fixFile = $<
                                                                                                               
if ("" == $fixFile)then
  set fixFile = analysis/defaultMoreFixFile.cc
endif
if (!(-e $fixFile)) then
  echo "$fixFile does not exist. Exiting."
  exit
endif

# copy the files to the output directory for the record:
cp $floatFile $dir
cp $parInputFile $dir

echo -n "remove each job's csh file when done with it? [y/n] "
set removeCshFile = $<

echo -n "remove each job's log file when done with it? [y/n] "
set removeLogFile = $<

set cshFile = $dir/fitChunksToy-${stc}.csh
set logFile = $dir/fitChunksToy-${stc}.log

echo "#! /usr/local/bin/tcsh -f" > $cshFile
echo "bbrroot -l -b <<EOF " >> $cshFile

echo .x globals/setup.cc >> $cshFile
echo fitChunksToy\(2047, $firstJob, $lastJob, $stc, $etc,\"$dir\/\"\, \"$fixFile\"\, \"$parInputFile\"\, \"$floatFile\"\) >>$cshFile
echo ".q" >> $cshFile
echo "EOF" >> $cshFile

chmod +x $cshFile

if ("$queue" != "") then
   echo "submit job to queue $queue"
   bsub -q $queue -o $logFile $cshFile 
else

echo "Running file $cshFile at local machine"
  $cshFile >& $logFile
endif

if ($removeCshFile == "y") then
    rm -f $cshFile  
endif 

if ($removeLogFile == "y") then
    rm -f logFile  
endif 


  
