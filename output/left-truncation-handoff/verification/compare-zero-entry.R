# Usage: Rscript verification/compare-zero-entry.R
# Set R_LIBS first if the updated SurvPT is installed in a local library.
library(SurvPT)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this script using Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
baseline <- readRDS(file.path(script_dir, "left-truncation-before.rds"))
if (!"entry" %in% names(formals(survPT))) stop("Install the updated SurvPT first.")
for (key in names(baseline$cases)) {
    args <- baseline$cases[[key]]
    old <- baseline$fits[[key]]
    for (variant in list(list(), list(entry = 0), list(entry = rep(0, nrow(args$X))))) {
        current <- do.call(survPT, c(args, variant))
        if (!identical(old, current)) stop("Changed legacy output for case: ", key)
    }
}
cat(length(baseline$cases), "pre-change fits x 3 entry variants: identical() for every complete result.\n")
