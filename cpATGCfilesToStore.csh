#!/bin/tcsh

# make multiple directories in FNAL store/user/ area 

cd /uscmst1b_scratch/lpc1/3DayLifetime/kalanand/

foreach dirname ( `ls -1 *.root` ) 
	echo "submitting: srmcp $dirname lnujj/MCFM_WW_aTGC/Grid2D_LambdaZ_KappaGamma/$dirname" 
	srmcp -2 -debug=true "file://localhost//uscmst1b_scratch/lpc1/3DayLifetime/kalanand/$dirname"  "srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/11/store/user/lnujj/MCFM_WW_aTGC/Grid2D_LambdaZ_KappaGamma/$dirname"
         
        

##foreach dirname ( `ls -1` ) 
##	echo "submitting: srmcp $dirname lnujj/aTGC_Madgraph/wz/$dirname" 
##	srmcp -2 -debug=true "file://localhost//uscmst1b_scratch/lpc1/3DayLifetime/kalanand/wzTGC-MadGraph/$dirname"  "srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/11/store/user/lnujj/aTGC_Madgraph/wz/$dirname"
##
         
    end
    echo ""



