options(stringsAsFactors = FALSE)

message("Patching Phase13.7A script for empty/unreadable headers...")

script_file <- "71_data_analysis_reinforcement/phase13_7_instrument_transparency/phase13_7A_discover_instrument_files.R"

stopifnot(file.exists(script_file))

lines <- readLines(script_file, warn = FALSE, encoding = "UTF-8")

start <- grep("^get_header <- function", lines)[1]
end <- grep("^has_any <- function", lines)[1] - 1

if (is.na(start) || is.na(end) || end <= start) {
  stop("Could not locate function block to patch.")
}

new_block <- c(
  "get_header <- function(path) {",
  "  con <- NULL",
  "  out <- NA_character_",
  "  tryCatch({",
  "    con <- open_connection(path)",
  "    tmp <- readLines(con, n = 1, warn = FALSE)",
  "    if (length(tmp) == 0 || is.na(tmp[1]) || !nzchar(tmp[1])) {",
  "      out <- NA_character_",
  "    } else {",
  "      out <- tmp[1]",
  "    }",
  "  }, error = function(e) {",
  "    out <<- NA_character_",
  "  }, finally = {",
  "    if (!is.null(con)) try(close(con), silent = TRUE)",
  "  })",
  "  out",
  "}",
  "",
  "guess_sep <- function(header_line, path) {",
  "  if (length(header_line) != 1 || is.na(header_line) || !nzchar(header_line)) return(\"\\t\")",
  "  count_fixed <- function(text, pattern) {",
  "    hit <- gregexpr(pattern, text, fixed = TRUE)[[1]]",
  "    if (length(hit) == 1 && hit[1] == -1) return(0L)",
  "    length(hit)",
  "  }",
  "  n_tab <- count_fixed(header_line, \"\\t\")",
  "  n_comma <- count_fixed(header_line, \",\")",
  "  if (grepl(\"\\\\.csv(\\\\.gz)?$\", path, ignore.case = TRUE)) return(\",\")",
  "  if (n_comma > n_tab) return(\",\")",
  "  \"\\t\"",
  "}",
  "",
  "split_header <- function(header_line, sep) {",
  "  if (length(header_line) != 1 || is.na(header_line) || !nzchar(header_line)) return(character(0))",
  "  out <- strsplit(header_line, split = sep, fixed = TRUE)[[1]]",
  "  out <- trimws(out)",
  "  out[nzchar(out)]",
  "}",
  ""
)

patched <- c(
  lines[seq_len(start - 1)],
  new_block,
  lines[(end + 1):length(lines)]
)

backup_file <- paste0(script_file, ".before_empty_header_patch")
file.copy(script_file, backup_file, overwrite = TRUE)

writeLines(patched, script_file, useBytes = TRUE)

message("Patch completed.")
message("Backup created: ", backup_file)
message("Patched script: ", script_file)
