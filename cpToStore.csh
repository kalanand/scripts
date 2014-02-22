#!/bin/tcsh

# copy data files from FNAL dCache resilient to FNAL store/user/ using srmcp

cd /pnfs/cms/WAX/resilient/kalanand/Madgraph_ttbarjets/

foreach datafile ( `ls *root` ) 
	echo "submitting: srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/resilient/kalanand/Madgraph_ttbarjets/$datafile  srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/11/store/user/kalanand/MC/1_6_12/ttbarjets/$datafile" 

	 srmcp -2 -debug=true "srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/resilient/kalanand/Madgraph_ttbarjets/$datafile" "srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/11/store/user/kalanand/MC/1_6_12/ttbarjets/$datafile"
    end
    echo ""



