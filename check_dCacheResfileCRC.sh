#!/bin/bash

if [ -z "$1" ] ; then 
  echo "please provide a path/file for checking" ; 
  exit -1 ; 
fi 


pathFile="`dirname $1`" ; 
lfn=`echo $1 | sed 's!^/resilient/!/!'` ; 
lfn=`echo $lfn | sed 's!^/pnfs.*/resilient/!/!'` ; 
# lfn="/resilient$lfn" ; 


echo $pathFile | grep "/pnfs/" >/dev/null 2>&1 ; 
if [ $? -eq 1 ] ; then pathFile="/pnfs/cms/WAX/resilient/`dirname $lfn`" ; fi ; 


lfnLocation=`curl -m 3600 -s http://cmspnfs1:1031/lfnLocationEcrcResilient.php?rpath=$lfn` ; 
cksumLen=`curl -m 3600 -s $lfnLocation 2>/dev/null | tail -1 | sed 's/.* (0x//' | sed 's/)//'` ; 
if [ -z "$cksumLen" ] ; then cksumLen=0 ; fi ; 
if [ ! -z "$2" ] ; then localFile="$2" ; else localFile=`basename $1` ; fi ; 


if [ ! -f "$pathFile/$localFile" ] ; then echo "Mistmacth: Local file not found" ; exit -2 ; fi ; 
cksumLocal=`cat "$pathFile/.(use)(2)($localFile)" 2>/dev/null | tail -1 | sed 's/.*c=//' | sed 's/;.*//' | sed 's/.*://'` ; 
cksumLocal2=`echo $cksumLocal | sed 's!^0*!!'`
if [ -z "$cksumLocal" ] ; then cksumLocal=0 ; fi ; 


if [ "$cksumLen" != "$cksumLocal" ] && [ "$cksumLen" != "$cksumLocal2" ] ; then
    echo "Mismatch: $cksumLen!=$cksumLocal or $cksumLocal2 " ;
    exit -1 ;
fi ;
echo "Match: $cksumLen=$cksumLocal " ; 
exit 0 ; 


