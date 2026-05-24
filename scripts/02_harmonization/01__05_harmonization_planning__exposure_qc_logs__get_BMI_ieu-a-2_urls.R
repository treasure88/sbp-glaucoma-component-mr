dataset_id <- "ieu-a-2"
log_dir <- "../../05_harmonization_planning/exposure_qc_logs"

if (!requireNamespace("ieugwasr", quietly = TRUE)) {
  install.packages("ieugwasr", repos = "https://cloud.r-project.org")
}

library(ieugwasr)

if (!nzchar(Sys.getenv("OPENGWAS_JWT"))) {
  stop("OPENGWAS_JWT is not set")
}

files <- ieugwasr::gwasinfo_files(dataset_id)
char_vals <- unlist(files[, vapply(files, is.character, logical(1)), drop = FALSE], use.names = FALSE)
char_vals <- char_vals[!is.na(char_vals)]

all_file <- file.path(log_dir, "BMI__ieu-a-2_all_opengwas_urls.txt")
vcf_file <- file.path(log_dir, "BMI__ieu-a-2_vcf_url.txt")
tbi_file <- file.path(log_dir, "BMI__ieu-a-2_tbi_url.txt")

writeLines(char_vals, all_file)

vcf_url <- char_vals[grepl("ieu-a-2\\.vcf\\.gz($|\\?)", char_vals)]
tbi_url <- char_vals[grepl("ieu-a-2\\.vcf\\.gz\\.tbi($|\\?)", char_vals)]

if (length(vcf_url) < 1) {
  stop("Could not find ieu-a-2.vcf.gz URL")
}

writeLines(vcf_url[1], vcf_file)

if (length(tbi_url) >= 1) {
  writeLines(tbi_url[1], tbi_file)
} else {
  writeLines("", tbi_file)
}

cat("Fetched fresh OpenGWAS temporary URLs for ", dataset_id, "\n", sep = "")
