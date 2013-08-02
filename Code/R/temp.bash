for file in *.R
do
	# (cat disclaimer.txt; cat $file) > $file.new
	mv $file.new $file
done