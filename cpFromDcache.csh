#!/bin/tcsh

# copy data files from FNAL dCache resilient to local area

cd /pnfs/cms/WAX/resilient/kalanand/ZeeJet_OctX_7TeV/

foreach datafile ( `ls *root` ) 
	echo "submitting: dccp $datafile  /uscms/home/kalanand/cms/jet/CMSSW_3_1_4/src/JetMETCorrections/ZJet/test/$datafile" 

	 dccp $datafile  /uscms/home/kalanand/cms/jet/CMSSW_3_1_4/src/JetMETCorrections/ZJet/test/$datafile
    end
    echo ""



