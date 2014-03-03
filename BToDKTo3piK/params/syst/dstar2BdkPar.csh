#! /usr/local/bin/tcsh -f
# Creates BDK .par files from the dstar .par file

foreach parFile ($argv)
    echo "Processing $parFile"
    set file = `basename $parFile .par`
    cat $parFile | grep -v norm | sed "s/dstar/dalitzHolderN.sigGoodD0/" > $file.Bdk.par
#    cat $parFile | grep -v norm | sed "s/dstar/dalitzHolderP.sigGoodD0/" >> $file.Bdk.par
end
