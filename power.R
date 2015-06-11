require(pwr, quietly=TRUE)

P.VALUE <- 0.05
ERROR <- 0.001
USAGE <- c(
  "Usage: $ Rscript power.R data.tsv",
  "       $ cat data.tsv | Rscript power.R",
  "Input file must be tab-separated: column 1 is the maf, column 2 is the coverage.",
  "The rest of the columns don't matter."
)

# Read arguments.
args <- commandArgs(trailingOnly=TRUE)
if ((length(args) >= 1 && args[1] == '-h') || length(args) > 1) {
  cat(paste(USAGE), sep="\n")
  quit()
}

if (length(args) == 1) {
  infile <- args[1]
} else {
  infile <- file('stdin')
}

data <- read.delim(infile, header=FALSE)
colnames(data) <- c('maf', 'cvg')

for (i in 1:nrow(data)) {
  maf <- data$maf[i]
  cvg <- data$cvg[i]
  power <- pwr.p.test(h=ES.h(maf, ERROR), n=cvg, sig.level=P.VALUE, alternative='greater')
  cat(power$power, sep="\n")
}
