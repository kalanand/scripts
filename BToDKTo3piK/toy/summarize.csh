#
# reads the .out files prduced by fitToy.cc and makes ntuples.
#

cd $1
rm -rf summarize.dat
echo -n " " > summarize.dat
@ n = 0
foreach file (*.out)
  cat $file |sed s/fitStatus// >> summarize.dat
  @ n = $n + 1
end

echo "wrote $n files to summarize.dat"

