# $Id: copy.csh,v 1.4 2005/10/11 18:57:36 abi Exp $
#
echo "Run this script from workdir to copy BToDKTo3piK subdirectories to workdir"
echo " "

  set pwd=`pwd`
  if ("${pwd:t}" != "BToDKTo3piK" && "${pwd:t}" != "workdir") then
    echo ""
    echo " Abort. You have to do be in BToDKTo3piK or workdir to run this script."
    echo ""
    exit
  endif


set list = "analysis globals params utils kumac syst toy"
set copyAll = 0
foreach dir ($list)
  if (0 == $copyAll) then
    echo -n "copy  ../BToDKTo3piK/$dir/ to workdir ? [y/n, a to copy all] "
    set answer = $<
  endif

  if ("a" == "$answer") then
    set copyAll = 1
  endif

  if ("y" == "$answer" || 1 == $copyAll) then
    echo copying $dir/
    cp -rf ../BToDKTo3piK/$dir ../workdir
  endif
end
