sed 's/\*//g' $1 | awk '{print $1" "$3" "$2}' > $1.sed
