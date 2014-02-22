#!/bin/tcsh

# copy data files from FNAL dCache to CERN castor using srmcp

cd /pnfs/cms/WAX/resilient/kalanand/Madgraph_ttbarjets/

foreach datafile ( `ls *root` ) 
	echo "submitting: srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/resilient/kalanand/Madgraph_ttbarjets/$datafile  srm://srm-cms.cern.ch:8443//srm/managerv2?SFN=/castor/cern.ch/user/k/kalanand/Madgraph_ttbarjets/$datafile" 

	 srmcp -2 -debug=true "srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/resilient/kalanand/Madgraph_ttbarjets/$datafile" "srm://srm-cms.cern.ch:8443//srm/managerv2?SFN=/castor/cern.ch/user/k/kalanand/Madgraph_ttbarjets/$datafile"
    end
    echo ""



