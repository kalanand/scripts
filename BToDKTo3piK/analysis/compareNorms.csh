#
# This script compares the normalization parameters of a dalitzHolderN
# and a dalitzHolderP PDF. Put the normalization parameters in the files
# n.par and p.par, then run this script in the same directory.
# It prints the N and P parameters followed by their "asymmetry"
# (P-N)/(P+N).
#

set length = `wc p.par |awk '{print $1}'`
set i = 0
while ($i < $length)
  @ i = $i + 1
  set nn = `tail -$i n.par |head -1 |awk '{print $3}'`
  set pp = `tail -$i p.par |head -1 |awk '{print $3}'`
  echo -n "$nn $pp \t"
  calc "($pp-($nn))/($pp+($nn))"
end

