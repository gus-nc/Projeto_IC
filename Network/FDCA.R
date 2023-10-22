###
### Perform a FDCA as inspired by Tinker et al 2012
###-----Set Directory-----------------------------------------------------------
setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Laptop
setwd("D:/Drive/Other computers/Meu laptop/Code/Projeto_IC") # PC
###-----Define Functions--------------------------------------------------------
### Function to import matrices in the appropriate format
# x and e must be the same length
import_and_cut = function(x, e, n) {
  if (length(x) != length(e)) {
   stop("Both inputs must be the same length.")
  }
  webs = list()
  rep = 0
  for (i in x) {
    rep = rep + 1
    webs[[envs[rep]]] = t(as.matrix(read.csv(i, header = T, row.names = 1)))
    
    # Filter the matrix
    cols_to_remove = which(colSums(webs[[envs[rep]]]) < n) # Only for freq bc it reveals the discrete number of preys
    # Min of 5 prey items p/ species
    if (length(cols_to_remove) > 0) { # Check if is a int(0)
      webs[[envs[rep]]] = webs[[envs[rep]]][, -cols_to_remove] # Remove columns
    }
    rows_to_remove = which(rowSums(webs[[envs[rep]]]) == 0) # Only for freq bc it reveals the discrete number of preys
    # Min of 5 prey items p/ species
    if (length(rows_to_remove) > 0) { # Check if is a int(0)
      webs[[envs[rep]]] = webs[[envs[rep]]][-rows_to_remove,] # Remove rows
    }
  }
  
  return(webs)
} # Import a matrix and cut by n items for a species

### Function for executing an FDCA

fraction_analysis = function(x, fraction) {
  mat_fdca = list()
  mat_names = names(x)
  for (f in fraction) {
  rep = 0
  for (k in x) {
    rep = rep + 1
    k = apply(k, MARGIN = 2, function(x) sort(x, decreasing = TRUE))
    for (j in 1:dim(k)[2]) {
      data = k[,j] # Extract a column to analyze
      total_sum = sum(data)  # Calculate the total sum of the column
      threshold = f * total_sum # Determine the threshold for retaining values
      sorted_data = sort(data, decreasing = TRUE) # Sort the data in descending order
       
      # Initialize variables to store retained values and cumulative sum
      retained_values = numeric(0)
      cumulative_sum = 0
      
      # Iterate through the sorted data to retain values
      for (value in sorted_data) {
        if (cumulative_sum + value <= threshold) {
          retained_values = c(retained_values, value)
          cumulative_sum = cumulative_sum + value
        } 
        else { # Iterate one time to include thye category tha surpass the treshold
          retained_values = c(retained_values, value) 
          cumulative_sum = cumulative_sum + value
          break
        }
      }
      # Store only the retained_values in the new matrix
      for (i in 1:length(retained_values)) {
        k[i,j] = retained_values[i] 
      }
      for (i in (length(retained_values)+1):dim(k)[1]) {
        k[i,j] = 0
      }
    }
    mat_fdca[[mat_names[rep]]][[as.character(f)]] = k
  }
  }
  # Return retained values and their sum
  return(mat_fdca)
}

### Function for calculating Nestedness NODF for 1 matrix
Nestedness_NODF = function(x) {
  
  NODFc_temp = 0 # Create a temporary object for NODFc
  for (i in 1:(ncol(x)-1)) { # for every col i:n-1
    for (j in (i+1):ncol(x)) { # and for every col i+1:n
      if(sum(x[,i]) <= sum(x[,j])) { # if F(ci) <= F(cj)
        N_paired = 0 # The Nested Value for this pair of cols
      }
      else {
        PO = 0 # number of 1's in ci = cj, if cj != 0?
        N_j = sum(x[,j] == 1) # Cells > 0 in cj
        for (r in 1:nrow(x)) { # for every row in t cols sum PO
          if (x[r,i] == 1 && x[r,j] == 1){
            PO = PO+1
          }
        }
        N_paired = PO/N_j # The Nested Value for this pair of cols
      }
      NODFc_temp = NODFc_temp + N_paired
    }
  }
  NODFc = NODFc_temp*100
  # NODFr values for rows
  NODFr_temp = 0 # Create a temporary object for NODFr
  for (i in 1:(nrow(x)-1)) { # for every col i:n-1
    for (j in (i+1):nrow(x)) { # and for every col i+1:n
      if(sum(x[i,]) <= sum(x[j,])) { # if F(ci) <= F(cj)
        N_paired = 0 # The Nested Value for this pair of cols
      }
      else {
        PO = 0 # number of 1's in ci = cj, if cj != 0?
        N_j = sum(x[j,] == 1)  #Cells > 0 in cj
        for (c in 1:ncol(x)) { # for every row in t cols sum PO
          if (x[i,c] == 1 && x[j,c] == 1){
            PO = PO+1
          }
        }
        N_paired = PO/N_j # The Nested Value for this pair of cols
      }
      NODFr_temp= NODFr_temp + N_paired
    }
  }
  NODFr = NODFr_temp*100
  # NODF
  m = dim(x)[1] # Number of rows
  n = dim(x)[2] # Number of cols
  NODF_obs = 2*(NODFc+NODFr)/ # The final average of NODF for the env matrix
    (m*(m-1)+(n*(n-1)))
  return(c(NODFc, NODFr, NODF_obs))
}

###-----Import Martrices--------------------------------------------------------
# Import the matrices and create a webs list
# A matrix as > higher trophic level species=anurans (columns) 
# and lower trophic level species=prey (rows).
files = c("Diet_E_Freq_Xincomp.csv", "Diet_M_Freq_Xincomp.csv")
envs = c("freq_eu", "freq_mt")

net = import_and_cut(files, envs, n = 5)

###-----Execute FDCA------------------------------------------------------------
frac = seq(0.1, 0.9, by = 0.1)
frac_matrices = fraction_analysis(net, fraction = frac)

# Re-order the interaction matrices, removes rows that sum == 0 and make then binary
for (env in envs) {
  for (f in frac) {
    frac_matrices[[env]][[as.character(f)]] = ifelse(frac_matrices[[env]][[as.character(f)]] > 0, 1, frac_matrices[[env]][[as.character(f)]])
    
    row_order = order(rowSums(frac_matrices[[env]][[as.character(f)]]), decreasing = TRUE) # Calculate row sums
    frac_matrices[[env]][[as.character(f)]] = frac_matrices[[env]][[as.character(f)]][row_order, ] # Use the index vector to reorder the rows of the matrix
    col_order = order(colSums(frac_matrices[[env]][[as.character(f)]]), decreasing = TRUE) # Calculate col sums
    frac_matrices[[env]][[as.character(f)]] = frac_matrices[[env]][[as.character(f)]][,col_order] # Use the index vector to reorder the cols of the matrix
    rows_to_remove = which(rowSums(frac_matrices[[env]][[as.character(f)]]) == 0) # Only for freq bc it reveals the discrete number of preys
      if (length(rows_to_remove) > 0) { # Check if is a int(0
        frac_matrices[[env]][[as.character(f)]] =frac_matrices[[env]][[as.character(f)]][-rows_to_remove,] # Remove rows
      }
  }
}

# Incorportate the Original matrix
frac = c(frac, 1)

for (env in envs) {
  frac_matrices[[env]][["1"]] = net[[env]]
  frac_matrices[[env]][["1"]] = ifelse(frac_matrices[[env]][["1"]] > 0, 1, frac_matrices[[env]][["1"]])
  
  row_order = order(rowSums(frac_matrices[[env]][["1"]]), decreasing = TRUE) # Calculate row sums
  frac_matrices[[env]][["1"]] = frac_matrices[[env]][["1"]][row_order, ] # Use the index vector to reorder the rows of the matrix
  col_order = order(colSums(frac_matrices[[env]][["1"]]), decreasing = TRUE) # Calculate col sums
  frac_matrices[[env]][["1"]] = frac_matrices[[env]][["1"]][,col_order] # Use the index vector to reorder the cols of the matrix
  rows_to_remove = which(rowSums(frac_matrices[[env]][["1"]]) == 0) # Only for freq bc it reveals the discrete number of preys
  if (length(rows_to_remove) > 0) { # Check if is a int(0
    frac_matrices[[env]][["1"]] =frac_matrices[[env]][["1"]][-rows_to_remove,] # Remove rows
  }
  
}


###----- Nestedness analysis----------------------------------------------------
NODF_results = list()

for (env in envs) {
  for (f in frac) {
    NODF_results[[env]][[as.character]] = Nestedness_NODF(frac_matrices[[env]][[as.character(f)]])
  }
}

