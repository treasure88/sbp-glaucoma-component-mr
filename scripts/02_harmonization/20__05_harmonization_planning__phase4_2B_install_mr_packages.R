options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("===== Phase 4.2B package installation =====\n")
cat("R version:", as.character(getRversion()), "\n\n")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing CRAN package:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat("Already installed:", pkg, "\n")
  }
}

install_if_missing("remotes")

cat("\n===== Installing MRPRESSO from GitHub =====\n")
if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  tryCatch({
    remotes::install_github("rondolab/MR-PRESSO", upgrade = "never", dependencies = TRUE)
  }, error = function(e) {
    cat("ERROR installing MRPRESSO:\n")
    cat(conditionMessage(e), "\n")
  })
} else {
  cat("MRPRESSO already installed\n")
}

cat("\n===== Installing mr.raps from GitHub =====\n")
if (!requireNamespace("mr.raps", quietly = TRUE)) {
  tryCatch({
    remotes::install_github("qingyuanzhao/mr.raps", upgrade = "never", dependencies = TRUE)
  }, error = function(e) {
    cat("ERROR installing mr.raps:\n")
    cat(conditionMessage(e), "\n")
  })
} else {
  cat("mr.raps already installed\n")
}

cat("\n===== Final package availability =====\n")
cat("MRPRESSO:", requireNamespace("MRPRESSO", quietly = TRUE), "\n")
cat("mr.raps:", requireNamespace("mr.raps", quietly = TRUE), "\n")
