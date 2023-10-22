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


###-----Import Martrices--------------------------------------------------------
# Import the matrices and create a webs list
# A matrix as > higher trophic level species=anurans (columns) 
# and lower trophic level species=prey (rows).
files = c("Diet_E_Freq_Xincomp.csv", "Diet_M_Freq_Xincomp.csv")
envs = c("freq_eu", "freq_mt")

net = import_and_cut(files, envs, n = 5)

###-----Execute FDCA------------------------------------------------------------
frac = c(0.1, 0.9)
frac_matrices = fraction_analysis(net, fraction = frac)


# Re-order the interaction matrices
for (env in envs) {
  row_order = order(rowSums(webs[[env]]), decreasing = TRUE) # Calculate row sums
  webs[[env]] = webs[[env]][row_order, ] # Use the index vector to reorder the rows of the matrix
  col_order = order(colSums(webs[[env]]), decreasing = TRUE) # Calculate col sums
  webs[[env]] = webs[[env]][,col_order] # Use the index vector to reorder the cols of the matrix
}
