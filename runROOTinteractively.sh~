#!/bin/tcsh

# convert all eps figures in this directory to pdf

foreach datafile ( `ls *.eps` ) 
        set fname = `echo $datafile | awk '{sub(".eps",".pdf"); print $0}'` 
	echo "submitting: convert $datafile $fname" 
        convert $datafile $fname	 
    end
    echo ""



