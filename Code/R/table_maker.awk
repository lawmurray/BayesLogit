BEGIN {
    printf "\\begin{tabular}{l r r r r r r r} \n"
    printf("%9s ", "Method")
}
{
    if (NR==1) for (i=1; i<=9; i++) printf(" & %16s", $i)
    else {
	printf " \\\\ \n"
	printf ("%9s ", $1)
	for (i=2; i<=10; i++) {
	    printf (" & %16.2f", $i) ;
	}
    }

}
END {
    print "\n \\end{tabular}"
}