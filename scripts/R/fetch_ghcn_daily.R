### download the data from NOAA at
### ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/
# !/usr/bin/env Rscript
library("optparse")
library("data.table")
library("tools")
library("R.utils")

fetch_file <- function(url, dest, overwrite = FALSE) {
  if (file.exists(dest)) {
    cat(dest, " exists already\n")
  }
  if (overwrite || (!file.exists(dest))) {
    download.file(url, dest)
  }
}
fselect <- function(source) {
  # select columns from .csv.gz file and creat .Rdata object
  t0 <- Sys.time()
  dt <- fread(source)
  dt <- setNames(dt, c("ID", "Date", "Element", "Value", "V5", "V6", "V7", "V8"))
  cat("Lines ", nrow(dt), "")
  t1 <- Sys.time()
  cat("t1 =", t1 - t0, "\n")
  # Select Only rows with TMIN and TMAX from dt$V3
  dt <- dt[which(dt$Element %in% c("TMIN", "TMAX")), ]
  cat("selected ", nrow(dt), "")
  t2 <- Sys.time()
  cat("t2 =", t2 - t1, "\n")
  # Create an array of [Date][Type][Element]
  # : chr [1:365] "20190101" "20190102" "20190103" "20190104" ...
  # : ID     : chr [1:14226] "AE000041196"
  # : Element: chr [1:2] "TM
  w <- tapply(dt$Value, dt[, c("Date", "ID", "Element")], as.integer)
  cat("dim =", dim(w), "")
  t3 <- Sys.time()
  cat("t3 =", t3 - t2, "\n")
  # TODO: Returns TMIN+TMAX? WHY does this make sense?
  return((w[, , 1] + w[, , 2]))
}

create_R_obj <- function(source, overwrite) {
  file_path_without_gz <- file_path_sans_ext(source)
  file_path_without_csv <- file_path_sans_ext(file_path_without_gz)
  dest <- paste0(file_path_without_csv, ".RData")
  if (file.exists(dest)) {
    cat(dest, " exists already\n")
  }
  if (overwrite | !file.exists(dest)) {
    cat("creating ", dest, "\n")
    filename <- basename(file_path_without_csv)
    ## year <- gsub("[^0-9.-]", "", filename)
    robj <- filename
    assign(
      robj,
      fselect(source),
      envir = .GlobalEnv
    )
    save(
      list = robj,
      file = dest,
      compress = "xz"
    )
    cat("done!\n")
  }
}

# Max time until download interrupted
options(timeout = 400)
ftp_path <- "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily"

option_list <- list(
  make_option(c("-p", "--ftp"),
    type = "character", default = ftp_path,
    help = "ftp server address [default= %default]", metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = "../../data/input/ghcn",
    help = "output folder name without year [default= %default]", metavar = "character"
  ),
  make_option(c("-y", "--year"),
    type = "integer", default = 2019,
    help = "Year [default= %default]", metavar = "number"
  ),
  make_option(c("-s", "--stations"),
    type = "logical", default = TRUE,
    help = "download stations ghcnd-stations.txt [default= %default]",
    metavar = "logical"
  ),
  make_option(c("-r", "--readme"),
    type = "logical", default = TRUE,
    help = "download readme.txt [default= %default]", metavar = "logical"
  ),
  make_option(c("c", "--convert"),
    type = "logical", default = TRUE,
    help = "Convert csv.gz to .Rdata file [default= %default]", metavar = "logical"
  ),
  make_option(c("--overwrite"),
    type = "logical", default = FALSE,
    help = "overwrite existing files [default= %default]", metavar = "logical"
  ),
  make_option(c("-v", "--verbose"),
    type = "logical", default = TRUE,
    help = "verbose output [default= %default]", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ================== START REAL PROGRAM =========================
if (sys.nframe() == 0) {
  # runs only when script is run by itself
  opt$output <- file.path(opt$output, opt$year)
  dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
  fetch_file(
    file.path(ftp_path, paste0("by_year/", opt$year, ".csv.gz")),
    file.path(opt$output, paste0("d", opt$year, ".csv.gz")),
    opt$overwrite
  )
  if (opt$readme) {
    fetch_file(
      file.path(ftp_path, "readme.txt"),
      file.path(opt$output, "readme.txt"),
      opt$overwrite
    )
  }
  if (opt$stations) {
    fetch_file(
      file.path(ftp_path, "ghcnd-stations.txt"),
      file.path(opt$output, "ghcnd-stations.txt"),
      opt$overwrite
    )
  }
  if (opt$convert) {
    create_R_obj(file.path(opt$output, paste0("d", opt$year, ".csv.gz")), opt$overwrite)
  }
}
