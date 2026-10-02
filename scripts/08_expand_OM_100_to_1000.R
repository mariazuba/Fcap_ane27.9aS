# ============================================================
# 08_expand_OM_100_to_1000.R
# Expand conditioned OMs from 100 to 1000 simulation trajectories
#
# Stock / analysis:
#   European anchovy, ane.27.9aS (Gulf of Cádiz)
#
# Role in the workflow:
#   This script follows the construction of the 100 conditioned
#   stock-recruitment reference objects in 07_build_SR_reference.R.
#   It expands each of the 100 historical OM conditions to 10
#   replicate trajectories for subsequent stochastic projection.
#
# Purpose:
#   - Load the 100-iteration FLBiol and FLSRsim reference objects.
#   - Build the explicit 100 OM x 10 replicate expansion map.
#   - Expand FLQuant and FLBiol objects to 1000 iterations.
#   - Expand the conditioned FLSRsim object consistently.
#   - Preserve the correspondence among trajectory, historical OM,
#     bootstrap condition, BH parameters and sigmaR.
#
# Inputs:
#   data/Rdata/biols_conditioned_100iter.rds
#   data/Rdata/SRs_reference_conditioned_100iter.rds
#   data/Rdata/SR_reference_iteration_map.csv
#
# Outputs:
#   data/Rdata/biols_conditioned_1000iter.rds
#   data/Rdata/SRs_reference_conditioned_1000iter.rds
#   data/Rdata/iteration_map_1000.csv
#
# Simulation structure:
#   Historical OMs:       100
#   Replicates per OM:     10
#   Total trajectories:  1000
#
# Important:
#   This script introduces no new stochastic variation.
#   The 10 replicates of each OM inherit the same conditioned
#   biological and stock-recruitment state.
#   Recruitment process error is introduced in a subsequent step.
# ============================================================

# ============================================================
# 0. PACKAGES
# ============================================================

rm(list = ls())

library(FLBEIA)
library(FLCore)
library(dplyr)
library(here)

# ============================================================
# 1. SETTINGS AND PATHS
# ============================================================

rds_boot_dir <- file.path(here(),"data","Rdata")

n_om <- 100
n_rep <- 10
n_it <- n_om * n_rep

# ============================================================
# 2. LOAD INPUT OBJECTS
# ============================================================
# Load the 100 conditioned biological OMs, their corresponding
# stock-recruitment objects and the validated OM-SR reference map.

biols_100 <- readRDS(file.path(rds_boot_dir,"biols_conditioned_100iter.rds"))
SRs_100 <- readRDS(file.path(rds_boot_dir,"SRs_reference_conditioned_100iter.rds"))
sr_reference_map <- read.csv(file.path(rds_boot_dir,"SR_reference_iteration_map.csv"))

# ============================================================
# 3. VALIDATE INPUT STRUCTURE
# ============================================================
# Confirm that the reference map contains one complete SR parameter
# set for each of the 100 historical OMs before expansion.

stopifnot(nrow(sr_reference_map) == n_om)
stopifnot(dplyr::n_distinct(sr_reference_map$iter) == n_om)
stopifnot(!anyNA(sr_reference_map$bootstrap))
stopifnot(!anyNA(sr_reference_map$a))
stopifnot(!anyNA(sr_reference_map$b))
stopifnot(!anyNA(sr_reference_map$sigmaR))

# ============================================================
# 4. BUILD 100 -> 1000 EXPANSION MAP
# ============================================================
# Assign 10 replicate trajectories to each historical OM and attach
# the corresponding bootstrap-specific SR parameters and sigmaR.

expansion_map <- data.frame(
  iter1000 = seq_len(n_it),
  om = rep(seq_len(n_om),each = n_rep),
  replicate = rep(seq_len(n_rep),times = n_om)
) %>%
  dplyr::left_join(sr_reference_map,by = c("om" = "iter"))

source_iter <- expansion_map$om

# ============================================================
# 5. VALIDATE EXPANSION MAP
# ============================================================
# Verify the intended 100 x 10 design and confirm that every expanded
# trajectory retains a complete link to its source OM and SR condition.

stopifnot(nrow(expansion_map) == n_it)
stopifnot(dplyr::n_distinct(expansion_map$iter1000) == n_it)
stopifnot(dplyr::n_distinct(expansion_map$om) == n_om)
stopifnot(all(table(expansion_map$om) == n_rep))
stopifnot(all(table(expansion_map$replicate) == n_om))
stopifnot(!anyNA(expansion_map$bootstrap))
stopifnot(!anyNA(expansion_map$a))
stopifnot(!anyNA(expansion_map$b))
stopifnot(!anyNA(expansion_map$sigmaR))

# ============================================================
# 6. HELPER FUNCTION — EXPAND FLQuant
# ============================================================
# Expand the iteration dimension of an FLQuant using source_iter.
# Single-iteration objects are replicated across all target iterations;
# multi-iteration objects inherit values from their source OM.

expand_FLQuant <- function(x,source_iter) {
  if(dim(x)[6] == 1) {
    out <- x[,,,,,rep(1,length(source_iter)),drop = FALSE]
  } else {
    out <- x[,,,,,source_iter,drop = FALSE]
  }

  dimnames(out)$iter <- as.character(seq_along(source_iter))

  out
}

# ============================================================
# 7. HELPER FUNCTION — EXPAND FLBiol
# ============================================================
# Expand FLQuant slots and the FLQuant components contained in
# predictModel slots while preserving the original FLBiol structure.

expand_FLBiol <- function(x,source_iter) {
  out <- x

  for(s in slotNames(x)) {
    z <- slot(x,s)

    if(inherits(z,"FLQuant")) {
      slot(out,s,check = FALSE) <- expand_FLQuant(z,source_iter)
    }

    if(inherits(z,"predictModel")) {
      z_new <- z
      z_new@.Data[[1]] <- expand_FLQuant(z@.Data[[1]],source_iter)
      slot(out,s,check = FALSE) <- z_new
    }
  }

  validObject(out)

  out
}

# ============================================================
# 8. EXPAND FLBiol TO 1000 TRAJECTORIES
# ============================================================
# Replicate each conditioned biological OM according to source_iter.

biols_1000 <- biols_100
biols_1000$ANE <- expand_FLBiol(biols_100$ANE,source_iter)

# ============================================================
# 9. EXPAND FLSRsim TO 1000 TRAJECTORIES
# ============================================================
# Expand the dynamic FLSRsim slots and BH parameter array using the
# same source-OM mapping applied to the biological object.

SRs_1000 <- SRs_100
SRs_1000$ANE <- SRs_100$ANE

SRs_1000$ANE@rec <- expand_FLQuant(SRs_100$ANE@rec,source_iter)
SRs_1000$ANE@ssb <- expand_FLQuant(SRs_100$ANE@ssb,source_iter)
SRs_1000$ANE@uncertainty <- expand_FLQuant(SRs_100$ANE@uncertainty,source_iter)
SRs_1000$ANE@proportion <- expand_FLQuant(SRs_100$ANE@proportion,source_iter)
SRs_1000$ANE@params <- SRs_100$ANE@params[,,,source_iter,drop = FALSE]

dimnames(SRs_1000$ANE@params)$iter <- as.character(seq_len(n_it))

# ============================================================
# 10. VALIDATE EXPANDED OBJECTS
# ============================================================
# Confirm that all expanded objects contain 1000 trajectories and
# that selected trajectories reproduce their intended source-OM
# stock-recruitment parameters exactly.

stopifnot(dim(biols_1000$ANE@n)[6] == n_it)
stopifnot(dim(SRs_1000$ANE@rec)[6] == n_it)
stopifnot(dim(SRs_1000$ANE@ssb)[6] == n_it)
stopifnot(dim(SRs_1000$ANE@uncertainty)[6] == n_it)
stopifnot(dim(SRs_1000$ANE@proportion)[6] == n_it)
stopifnot(dim(SRs_1000$ANE@params)[4] == n_it)

stopifnot(isTRUE(all.equal(SRs_1000$ANE@params[,,,1],SRs_100$ANE@params[,,,1])))
stopifnot(isTRUE(all.equal(SRs_1000$ANE@params[,,,10],SRs_100$ANE@params[,,,1])))
stopifnot(isTRUE(all.equal(SRs_1000$ANE@params[,,,11],SRs_100$ANE@params[,,,2])))
stopifnot(isTRUE(all.equal(SRs_1000$ANE@params[,,,1000],SRs_100$ANE@params[,,,100])))

# ============================================================
# 11. SAVE OUTPUT OBJECTS
# ============================================================
# Save the expanded biological and stock-recruitment objects together
# with the explicit 1000-trajectory mapping used by subsequent scripts.

saveRDS(biols_1000,file.path(rds_boot_dir,"biols_conditioned_1000iter.rds"))
saveRDS(SRs_1000,file.path(rds_boot_dir,"SRs_reference_conditioned_1000iter.rds"))
write.csv(expansion_map,file.path(rds_boot_dir,"iteration_map_1000.csv"),row.names = FALSE)
