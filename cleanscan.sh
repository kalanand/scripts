sed 's/\*//g' $1 | egrep -v R | egrep -v n | sed '/^$/d' | awk '!_[$4]++' | awk '{print $2" "$3" "$4" "$5" "$6}' > $1.txt
