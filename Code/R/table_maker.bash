for table in "$@"
do
	awk -f table_maker.awk $table
done
