#!/usr/bin/env bash
# $Id: normalizeAll.sh,v 1.1 2008/01/02 20:45:55 fwinkl Exp $
# Merge the normalization constants written by normalizeAll.cc with
# original par files

list="ksppPdf kkpPdf kskkPdf kskpPdf"
normPat="normRe|normIm|normD|x0|y0"

for f in $list; do
    cat "../BToDKTo3piK/params/$f.par" | egrep -v $normPat > "$f.par"
    cat "$f" | egrep $normPat >> "$f.par"
done
