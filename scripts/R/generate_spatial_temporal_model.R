# generate spatial-temporal matrices
# & data files
library("INLA")
library("optparse")
library("zeallot")
# source files with write functions
source("./fun_write_file.R")

generate_mesh <- function(spatial_res_param, ...) {
  return(inla.mesh.create(globe = spatial_res_param, ...))
}
create_spatial_mat <- function(mesh, base_path, order = 2) {
  # model order 2 for spatial
  sfem <- inla.mesh.fem(mesh, order = order)

  c0 <- as(sfem$c0, "dgCMatrix")
  mat_to_file_sym.fun(c0, "c0", base_path)
  g1 <- as(sfem$g1, "dgCMatrix")
  mat_to_file_sym.fun(g1, "g1", base_path)
  g2 <- as(sfem$g2, "dgCMatrix")
  mat_to_file_sym.fun(g2, "g2", base_path)

  return(sfem)
}
create_spatio_temporal_mat <- function(mesh, nt, base_path, order = 3) {
  # model order 3 for spatial-temporal
  ### define the spatial Finite Element Matrices
  sfem <- inla.mesh.fem(gmesh, order = order)
  ### define the temporal mesh for a given number of times
  tmesh <- inla.mesh.1d(loc = 1:nt)

  ### define the temporal Finite Element Matrices
  tfem <- inla.mesh.fem(tmesh, order = 2)

  # spatial matrices
  c0 <- as(sfem$c0, "dgCMatrix")
  mat_to_file_sym.fun(c0, "c0", base_path)
  g1 <- as(sfem$g1, "dgCMatrix")
  mat_to_file_sym.fun(g1, "g1", base_path)
  g2 <- as(sfem$g2, "dgCMatrix")
  mat_to_file_sym.fun(g2, "g2", base_path)
  g3 <- as(sfem$g3, "dgCMatrix")
  mat_to_file_sym.fun(g3, "g3", base_path)

  # temporal matrices
  # convert to dgCMatrix format
  M0_ddi <- tfem$c0
  # M0 <- sparseMatrix(i=c(1:M0_ddi@Dim[1]),j=c(1:M0_ddi@Dim[2]),x=M0_ddi@x,dims=list(M0_ddi@Dim[1],M0_ddi@Dim[2]))
  M0 <- as(tfem$c0, "dgCMatrix")
  mat_to_file_sym.fun(M0, "M0", base_path)
  M1 <- sparseMatrix(
    i = c(1, nt),
    j = c(1, nt),
    x = c(0.5, 0.5)
  )
  mat_to_file_sym.fun(M1, "M1", base_path)
  M2 <- tfem$g1
  mat_to_file_sym.fun(M2, "M2", base_path)
  return(list(sfem, tmesh, tfem))
}

##############################################################################
###########################     Main Program     #############################
##############################################################################
option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = "../../data/input/ghcn/2019",
    help = "Input file path [default= %default]", metavar = "character"
  ),
  make_option(c("-f", "--file"),
    type = "character", default = "d2019.RData",
    help = "Input file name [default= %default]", metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = "same as --input",
    help = "output folder directory [default= %default]", metavar = "character"
  ),
  make_option(c("-w", "--write"),
    type = "logical", default = TRUE,
    help = "Write to file [default= %default]", metavar = "logical"
  ),
  make_option(c("-s", "--spatial"),
    type = "integer", default = 2,
    help = "spatial resolution  ns=concat(spatial^2,2) [default= %default]", metavar = "character"
  ),
  make_option(c("-0", "--start"),
    type = "integer", default = 1,
    help = "start date [day], [default= %default].  1 = 1st of January", metavar = "number"
  ),
  make_option(c("-t", "--temporal"),
    type = "integer", default = 3,
    help = "temporal resolution [days], [nt < 366]. if t=1 => spatial model [default= %default]", metavar = "number"
  ),
  make_option(c("-v", "--verbose"),
    type = "logical", default = FALSE,
    help = "verbose output [default= %default]", metavar = "character"
  ),
  make_option(c("--solve"),
    type = "logical", default = FALSE,
    help = "solve system [default= %default]", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if (opt$output == "same as --input") {
  opt$output <- opt$input
}
#############################################################
nt <- opt$temporal
if (nt == 1) {
  SPATIAL_MODEL <- TRUE
  SPATIO_TEMPORAL_MODEL <- FALSE
  opt$output <- file.path(opt$input, "spatial")
  dir.create(opt$output)
} else {
  opt$output <- file.path(opt$input, "spatio_temporal")
  dir.create(opt$output)
  SPATIAL_MODEL <- FALSE
  SPATIO_TEMPORAL_MODEL <- TRUE
}
############################################################
t_resolution_list <- (opt$start):(opt$start + nt - 1)

################### MESH PART ##############################
### create a mesh over the globe for a given resolution
gmesh <- generate_mesh(opt$spatial)
### number of nodes in the spatial mesh
(ns <- gmesh$n)

########################### Create Mesh Matrices #############################

if (opt$write == TRUE && SPATIAL_MODEL) {
  dir_name <- file.path(opt$output, "spatial", paste0("ns", toString(ns)))
  dir.create(dir_name, recursive = TRUE)
  sfem <- create_spatial_mat(gmesh, dir_name)
} else {
  dir_name <- file.path(opt$output, paste0("ns", toString(ns), "_nt", toString(nt)))
  dir.create(dir_name, recursive = TRUE)
  c(sfem, tmesh, tfem) %<-% create_spatio_temporal_mat(gmesh, nt, dir_name)
}

########################### DATA PART ##########################################
# load data which in our case is the d19.RData
filepath <- file.path(opt$input, opt$file)
assign("data", get(load(filepath)))

# wo welche zeichen eingelesen werden oder so.
# Lagged differences
# 11  9 10  7 34  4  4  6
(ws <- diff(c(0, 11, 20, 30, 37, 71, 75, 79, 85)))

# read station information
# extract particular stations, longitute, latitude & altitude information
stations <- read.fwf(file.path(opt$input, "ghcnd-stations.txt"), ws[1:4])
colnames(stations) <- c("station", "latitude", "longitude", "altitude")
id.stations <- pmatch(colnames(data), stations$station)
if (opt$verbose) {
  head(stations)
  dim(stations)
  summary(stations)
  str(id.stations)
}
longlat <- as.matrix(stations[id.stations, 3:2])
### map the data locations (in longlat) to the sphere
coo.mollweide <- inla.mesh.map(longlat, "longlat", inverse = TRUE)
###
if (SPATIAL_MODEL) {
  dim(Ast <- inla.spde.make.A(gmesh, coo.mollweide))
} else {
  dim(Ast <- inla.spde.make.A(gmesh, kronecker(
    matrix(1, length(t_resolution_list)),
    coo.mollweide
  ),
  group = rep(t_resolution_list, each = nrow(coo.mollweide)),
  group.mesh = tmesh
  ))
}

if (opt$verbose) {
  summary(rowSums(Ast))
  summary(colSums(Ast))
}
###
### the observation vector /20 black magic in order to store in single precission
y <- as.vector(t(data[t_resolution_list, ])) / 20
table(y < (-90))
y[y < (-90)] <- NA

if (opt$verbose) {
  summary(y)
}
table(is.na(y))

table(y < (-50))
table(y > (50))

hist(y, -74:1531, xlim = c(-50, 50))

y[y < (-50)] <- NA
y[y > (50)] <- NA

length(y) == nrow(Ast)

length(iisel <- which(!is.na(y)))
A.select <- Ast[iisel, ]
y.select <- y[iisel]

if (opt$verbose) {
  summary(rowSums(A.select))
  summary(colSums(A.select))
}

### y = B * b + A * u + e
alt.km <- stations$altitude[id.stations] / 1000
summary(alt.km)
table(alt.km < (-0.95))
alt.km[alt.km < (-0.95)] <- 0 ### just for now
summary(alt.km)

# construct B
B <- cbind(
  1, ## intercept
  rep(alt.km, length(t_resolution_list))
)[iisel, ]
nb <- dim(B)[2]

# combine A.select & B
A.x <- cbind(A.select, B)

# store A.x and not A.st & B separately
if (opt$write) {
  mat_to_file_sp.fun(A.x, "Ax", opt$output)
  mat_to_file.fun(y.select, "y", opt$output)
}

# get number of observations
no <- length(y.select)


# ------------------------------------------------------------------------------------------- #
# for comparison
# need to make sure theta is the same !

### function to build Q.k(g.3) spatial matrix
Qgk.fun <- function(fem, g = 1, order = 3) {
  if (order == 1) {
    return(g^2 * fem$c0 + fem$g1)
  }
  if (order == 2) {
    return(g^4 * fem$c0 + 2 * g^2 * fem$g1 + fem$g2)
  }
  if (order == 3) {
    return(g^6 * fem$c0 + 3 * g^4 * fem$g1 +
      3 * g^2 * fem$g2 + fem$g3)
  }
  return(NULL)
}

Qst.fun <- function(sfem, tfem, g = exp(rep(0, 3))) {
  q1s <- Qgk.fun(sfem, g[3], 1)
  q2s <- Qgk.fun(sfem, g[3], 2)
  q3s <- Qgk.fun(sfem, g[3], 3)

  M0 <- tfem$c0
  nt <- nrow(M0)
  M1 <- sparseMatrix(
    i = c(1, nt),
    j = c(1, nt),
    x = c(0.5, 0.5)
  )
  M2 <- tfem$g1


  return(g[1] * (kronecker(M0, q3s) +
    kronecker(M1 * 2 * g[2], q2s) +
    kronecker(M2 * g[2]^2, q1s)))
}

if (TRUE) {
  theta <- c(5, -10, 2.5, 1)
  theta.u <- theta[2:4]

  Q.u <- Qst.fun(sfem, tfem, exp(theta.u))
  # mat_to_file_sym_CSR.fun(Q.u, "Qu_R", opt$output)
  # print("Q.u : ")
  # print(Q.u[1:10,1:10])


  Q.b <- sparseMatrix(i = c(1:nb), j = c(1:nb), x = rep(1e-5, nb), dims = list(nb, nb))

  Qub0 <- sparseMatrix(
    i = NULL,
    j = NULL,
    dims = c(nb, ns * nt)
  )

  Q.x <- rbind(
    cbind(Q.u, t(Qub0)),
    cbind(Qub0, Q.b)
  )

  Q.e <- Diagonal(length(y.select), exp(theta[1]))
  Q.xy <- Q.x + crossprod(A.x, Q.e) %*% A.x # crossprod = t(A)*Q.e (faster)

  # print(Q.xy[1:10,1:10])

  # L.xy <- t(chol(Q.xy))


  B.xey <- crossprod(A.x, Q.e) %*% y.select
  dim(B.xey)

  if (opt$write) {
    mat_to_file_sym.fun(Q.xy, "Qxy_R", opt$output)
    mat_to_file.fun(drop(B.xey), "bxy_R", opt$output)
  }

  if (opt$solve) {
    x <- solve(Q.xy, B.xey)
    mat_to_file.fun(drop(x), "x_sol_R", opt$output)
  }
}
