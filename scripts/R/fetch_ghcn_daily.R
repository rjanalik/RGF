### download the data from NOAA at
### ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/
# !/usr/bin/env Rscript
library("optparse")

fetch_file <- function(url, destintion, overwrite = FALSE) {
  if (file.exists(destintion)) {
    cat(destintion, " exists already\n")
  }
  if (overwrite | !file.exists(destintion)) {
    download.file(url, destintion)
  }
}
# Max time until download interrupted
options(timeout = 400)
ftp_path <- "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/"

option_list <- list(
  make_option(c("-p", "--ftp"),
    type = "character", default = ftp_path,
    help = "ftp server address [default= %default]", metavar = "character"
  ),
  make_option(c("-o", "--out"),
    type = "character", default = "../../data/output/",
    help = "output file name [default= %default]", metavar = "character"
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
  make_option(c("--overwrite"),
    type = "logical", default = FALSE,
    help = "overwrite existing files [default= %default]", metavar = "logical"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ================== START REAL PROGRAM =========================
fetch_file(
  paste0(ftp_path, "by_year/", opt$year, ".csv.gz"),
  paste0(opt$out, "d", opt$year, ".csv.gz"),
  opt$overwrite
)
if (opt$readme) {
  fetch_file(
    paste0(ftp_path, "readme.txt"),
    paste0(opt$out, "readme.txt"),
    opt$overwrite
  )
}
if (opt$stations) {
  fetch_file(
    paste0(ftp_path, "ghcnd-stations.txt"),
    paste0(opt$out, "ghcnd-stations.txt"),
    opt$overwrite
  )
}
