# --------- DIFFERENT FUNCTIONS TO WRITE TO FILE --------------------#
library(data.table)
options(scipen=999)

# write DENSE matrix to file
mat_to_file.fun <- function(M, file_name, file_path, append = FALSE, ns=NULL, nt=NULL){
  
  dim_M <- dim(M)

    if(is.null(dim_M)){
    dim_M <- c(length(M), 1)
  }
  
  M.df <- data.frame(M)
  
  if( is.null(ns) != TRUE & is.null(nt) != TRUE){
    input_file = file.path(file_path, paste(file_name, "_", toString(dim_M[1]), "_", toString(dim_M[2]), "_ns", toString(ns), "_nt", toString(nt), ".dat", sep=""))
  } else{
    input_file = file.path(file_path, paste(file_name, "_", toString(dim_M[1]), "_", toString(dim_M[2]), ".dat", sep=""))    
  }

  fwrite(M.df, input_file, append = append, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
  print(paste("wrote matrix to file under :", input_file))

}

# matrix has to be in sparse (dgC format), assumed to be quadratic
mat_to_file_sp.fun <- function(M, file_name, file_path, ns=NULL, nt=NULL){
  
  # drop numerical zero entries
  # ASSUMING DIAGONAL IS NONZERO!
  M <- drop0(M)

  dim_M <- M@Dim
  nnz <- nnzero(M) #+ sum(M@x==0) # have to also count diagonal entries
  M.df <- data.frame(c(dim_M[1], dim_M[2], nnz, M@i, M@p, M@x))
  
  if(is.null(ns) != TRUE & is.null(nt) != TRUE){
    output_file = file.path(file_path, paste(file_name, "_", toString(dim_M[1]), "_", toString(dim_M[2]), "_ns", toString(ns), "_nt", toString(nt), ".dat", sep=""))
  } else{
    output_file = file.path(file_path, paste(file_name, "_", toString(dim_M[1]), "_", toString(dim_M[2]), ".dat", sep=""))
  }

  fwrite(M.df, output_file, append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
  print(paste("wrote matrix to file under :", output_file))
  
}

# matrix has to be in sparse (dgC format) & symmetric
# pardiso wants upper triangular CSR format which is equivalent to
# lower triangular CSC format (which is this)
mat_to_file_sym.fun <- function(M, file_name, file_path){
  
  # drop numerical zero entries
  # ASSUMING DIAGONAL IS NONZERO!
  M <- drop0(M)

  M_lower <- tril(M)
  n <- M_lower@Dim[1]
  nnz <- nnzero(M_lower) #+ sum(M@x==0) # have to also count diagonal entries
  M_lower.df <- data.frame(c(n, n, nnz, M_lower@i, M_lower@p, M_lower@x))
  
  input_file = file.path(file_path, paste(file_name, "_", toString(n), ".dat", sep=""))
  fwrite(M_lower.df, input_file, append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
  print(paste("wrote matrix to file under :", input_file))
  
  return(n)
}

# write other constants to file
const_to_file.fun <- function(const, file_name, file_path){
  input_file = file.path(file_path, file_name)
  const.df <- data.frame((const))
  
  fwrite(const.df, input_file, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
  print(paste("wrote constants to file under :", input_file))
}
