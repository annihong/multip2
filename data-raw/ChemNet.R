
## code to prepare `ChemNet` dataset goes here
## code to prepare `ChemNet` dataset goes here
# influence reputation data
# square matrix with mutual influence attribution
# 1 = influential; 0 = not influential
# cells contain the ratings of row actors about column actors
infrep <- read.table(file=system.file("data-raw", "infrep.csv", package="multiplexP2"), header=T, row.names="label", sep=";")

# political/strategic information exchange data
# directed network
pol <- read.table(file=system.file("data-raw", "pol.csv", package="multiplexP2"), header=T, row.names="label", sep=";")

# scientific sender matrix
# row actor sends scientific/technical information to column actor
scifrom <- read.table(file=system.file("data-raw", "scifrom.csv", package="multiplexP2"), header=T, row.names="label", sep=";")

# scientific receiver matrix
# row actor receives scientific/technical information from column actor
scito <- read.table(file=system.file("data-raw", "scito.csv", package="multiplexP2"), header=T, row.names="label", sep=";")

intpos <- read.table(file=system.file("data-raw", "intpos.csv", package="multiplexP2"), header=T, row.names="label", sep=";")
# type of organization
types <- read.table(file=system.file("data-raw", "orgtypes.csv", package="multiplexP2"), header=T, row.names="label", sep=";")

# apply some changes to the data to make them network-compatible
govt <- types == "gov"
prefsim <- dist(intpos, method="euclidean", diag=F, upper=T) # equation 2
prefsim <- max(prefsim) - prefsim # equation 3
prefsim <- as.matrix(prefsim)
sci <- as.matrix(scito) * t(as.matrix(scifrom)) # equation 1 in the paper
pol <- as.matrix(pol) # convert to matrix object
infrep <- as.matrix(infrep) # convert to matrix object

ChemNet <- list(
  "infrep" = infrep,
  "pol" = pol,
  "sci" = sci,
  "prefsim" = prefsim,
  "govt" = govt
)
usethis::use_data(ChemNet, overwrite = TRUE)
