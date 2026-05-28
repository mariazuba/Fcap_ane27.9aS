# ============================================================
# Local runner for Fcap_ane27.9aS
# ============================================================

# Project root: use current working directory
root <- getwd()

# Source scripts/functions
source(file.path(root, "scripts", "data.R"))
source(file.path(root, "functions", "perfectObs4seas.R"))
source(file.path(root, "functions", "FcapBpaHCR_ane9aS.R"))
source(file.path(root, "scripts", "model.R"))
source(file.path(root, "scripts", "output_Fcap.R"))
source(file.path(root, "scripts", "output.R"))
source(file.path(root, "scripts", "report.R"))