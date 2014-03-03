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

echo -n "Which queue do you want to submit the job: "
set queue = $<

if ("" == $firstJob || "" == $lastJob) then
  echo empty argument given for either firstJob, lastJob. Exiting.
  exit
endif

echo -n "Enter parInputFile (path from workdir), <return> for default: "
set parInputFile = $<

if ("" == $parInputFile)then
  set parInputFile = toy/defaultParInputFile.cc
endif
if (!(-e $parInputFile)) then
  echo "$parInputFile does not exist. Exiting."
  exit
endif

echo -n "Enter floatFile (path from workdir), <return> for default: "
set floatFile = $<

if ("" == $floatFile)then
  set floatFile = toy/defaultFloatFile.cc
endif
if (!(-e $floatFile)) then
  echo "$floatFile does not exist. Exiting."
  exit
endif

# copy the files to the output directory for the record:
cp $floatFile $dir
cp $parInputFile $dir

echo -n "remove N1 from the fit? [y/n] "
set noN1 = $<

echo -n "remove N2 from the fit? [y/n] "
set noN2 = $<

echo -n "remove each job's csh file when done with it? [y/n] "
set removeCshFile = $<

echo -n "remove each job's log file when done with it? [y/n] "
set removeLogFile = $<


set job = $firstJob
while ($job <= $lastJob)
  set cshFile = $dir/fitToy-${job}.csh
  set logFile = $dir/fitToy-${job}.log
  rm -f $logFile

  rm -f cshFile
  echo "#! /usr/local/bin/tcsh -f" > $cshFile
  echo "bbrroot -l -b <<EOF " >> $cshFile

  echo .x globals/setup.cc >> $cshFile
  echo "doPlot = kFALSE" >> $cshFile

  if ($noN1 == "y") then
    echo "pdfOnResDK.useN1(kFALSE)"  >>$cshFile
  endif
  if ($noN2 == "y") then
    echo "pdfOnResDK.useN2(kFALSE)"  >>$cshFile
  endif


  echo fitToy\(\"${dir}/fitToy-${job}.out\",$job, \"$parInputFile\"\, \"$floatFile\"\) >>$cshFile

  echo ".q" >> $cshFile
  echo "EOF" >> $cshFile

  
  chmod +x $cshFile

  if ("$queue" != "") then
    echo "submit job to queue $queue"
    bsub -q $queue -o $dir/fitToy-${job}.log $cshFile 
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

  @ job = $job + 1
end

  
