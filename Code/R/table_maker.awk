## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

BEGIN {
    print  ARGV[1]
    printf "\\begin{tabular}{l r r r r r r r r } \n"
    printf("%16s ", "Method")
}
{
    if (NR==1) {
	for (i=1; i<=NF; i++) {
	    namewidth[i] = length($i) < 8 ? 8 : length($i) 
	    format = " & %" namewidth[i] "s"
	    sub(/\"/g, "", $i)
	    printf(format, $i)
	    ## printf(" & %16s", $i)
	}
    }
    else {
	printf " \\\\ \n"
	printf ("%16s ", $1)
	for (i=2; i<=NF; i++) {
	    format = " & %" namewidth[i-1] ".2f"
	    printf (format, $i) 
	    ## printf (" & %16.2f", $i)
	}
    }

}
END {
    print "\n \\end{tabular}"
}