#!/bin/tcsh

set DIR = "00_09 10_19 20_29 30_39 40_49 50_59 60_69 70_79 80_89 90_99"

    set MB = "0 1 2 3 4 5 6 7 8 9"
set FLIST = "00000 00001 00002 00003 00004 00005 00006 00007 00008 00009 00010 00011 00012 00013 00014 00015 00016 00017 0001800019 00020 00021 00022 00023 00024 00025 00026 00027 00028 00029 00030 00031 00032 00033 00034 00035 00036 00037 00038 00039 00040 00041 00042 00043 00044 00045 00046 00047 00048 00049 00050 00051 00052 00053 00054 00055 00056 00057 00058 00059 00060 00061 00062 00063 00064 00065 00066 00067 00068 00069 00070 00071 00072 00073 00074 00075 00076 00077 00078 00079 00080 00081 00082 00083 00084 00085 00086 00087 00088 00089 00090 00091 00092 00093 00094 00095 00096 00097 00098 0009"

#Let's copy all files from each of the above directories
foreach var1 ($DIR)
     foreach var2 ($MB)
         foreach var3 ($FLIST)
	     echo "submitting: srmcp file:////afs/cern.ch/user/k/kalanand/cms/softQCD_13TeV/float/$var1/pythia8_mb$var2_$var3.root  srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/11/store/user/lnujj/BOOST2012/softQCD_13TeV/float/$var1/pythia8_mb$var2_$var3.root" 
	     srmcp -2 -debug=true "file:////afs/cern.ch/user/k/kalanand/cms/softQCD_13TeV/float/$var1/pythia8_mb$var2_$var3.root" "srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/11/store/user/lnujj/BOOST2012/softQCD_13TeV/float/$var1/pythia8_mb$var2_$var3.root"
         end
     end
end

