cat("CESGA test OK\n")
cat("Working directory:\n")
print(getwd())

cat("\nR version:\n")
print(R.version.string)

cat("\nFiles in project:\n")
print(list.files())

cat("\nChecking key folders:\n")
print(dir.exists(c("scripts", "functions", "boot", "data", "outputs")))

cat("\nChecking input files:\n")
print(file.exists(c(
  "boot/data/stk_ane9aS.rds",
  "boot/data/ss3_ane9aS.rds",
  "scripts/data.R",
  "functions/perfectObs4seas.R",
  "functions/FcapBpaHCR_ane9aS.R"
)))