suppressPackageStartupMessages(library(INLA))
#inla.setOption(inla.timeout = 10*60)
inla.setOption(inla.timeout = 3*60*60)
source("../scripts/SPOTS/helpers.R")
source("../scripts/Highplex/helpers.R")
source("../INLA/LCAR.R")
source("../INLA/MCAR.R")
source("../INLA/CCAR.R")
source("../INLA/MCCAR.R")
source("../INLA/spotCCAR.R")
source("../INLA/spotMCCAR.R")

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 