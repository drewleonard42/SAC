for file in $@;
do 
	mv $file $file"_tmp"
done
for file in $@;
do 
	echo "Moving "$file
	suffix=${file#*.}
	prefix=${file%.*}
	prefix=$((10#$prefix + 1))
	newprefix=$(printf "%05i" $prefix)
	mv $file"_tmp" $newprefix"."$suffix
done
