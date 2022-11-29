#args[1] - name of the input file
#args[2] - name of the output file

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

in_table <- read.table (args [1], header=TRUE, row.names=1)

len <- dim (in_table) [1]

for (j in 1:len) {
    if (in_table [j, "Min"] < 0) {
        cat ("Error: invalid values specified\n", file=args[2], append=TRUE);
        stop ("ERROR: negative value of a minimum value for a parameter\n")
    }
    if (in_table [j, "Min"] > in_table [j, "Max"]) {
        cat ("Error: invalid values specified\n", file=args[2], append=TRUE);
        stop ("ERROR: minimum value exceeds maximum value for a parameter\n")
    }
}


for (j in 1:len) {
    if (in_table [j, "Min"] == in_table [j, "Max"]) {
        cat (row.names (in_table) [j], "\tfixed\t", in_table [j, "Min"], "\t1.0\n", sep="", file=args[2], append=TRUE)
    } else {
        gmean <- sqrt (in_table [j, "Min"] * in_table [j, "Max"])
        span <- log10 (in_table [j, "Max"] / gmean)
        cat (row.names (in_table) [j], "\tnorm\t", gmean, "\t", span, "\n", sep="", file=args[2], append=TRUE)
    }
}


